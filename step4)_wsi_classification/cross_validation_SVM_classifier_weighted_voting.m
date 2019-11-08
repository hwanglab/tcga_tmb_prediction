% this function is used for 5-folds cross valition by SVM classifier with linear kernel
% we test it for 3 class classifications g6, g7 and >=g8
% author: Hongming Xu, Cleveland Clinic, May 2018
% you may modify or use it, but you need give the credit for original
% author


function cross_validation_SVM_classifier_weighted_voting

clear vars;
% featPath={'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\2)mutation_prediction\low\','E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\2)mutation_prediction\mid\',...
%     'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\2)mutation_prediction\high\'};

featPath={'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\7)mutation_prediction_clustering\low\','E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\7)mutation_prediction_clustering\mid\',...
     'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\7)mutation_prediction_clustering\high\'};
 
Ftrain=[];
Ltrain=[];
Lwsi=[];
Lse=[];
Ww=[];

% step 1) feature organization
for ip=1:length(featPath)
    fp=featPath{ip};
    mats=dir(fullfile(fp,'*.mat'));
    
    for im=1:numel(mats)
        matfilename=strcat(fp,mats(im).name);
        load(matfilename);
        
        Ftrain=[Ftrain;feat_out(:,41:112)];
        Ww=[Ww;feat_out(:,end)];
        
        if ip==1
            Llow=ones(size(feat_out,1),1);
            Ltrain=[Ltrain;Llow];
            Lwsi=[Lwsi;1];
            Lse=[Lse;size(Ltrain,1)-size(feat_out,1)+1,size(Ltrain,1)];
        elseif ip==2
            Lmid=2*ones(size(feat_out,1),1);
            Ltrain=[Ltrain;Lmid];
            Lwsi=[Lwsi;2];
            Lse=[Lse;size(Ltrain,1)-size(feat_out,1)+1,size(Ltrain,1)];
        else
            Lhigh=3*ones(size(feat_out,1),1);
            Ltrain=[Ltrain;Lhigh];
            Lwsi=[Lwsi;3];
            Lse=[Lse;size(Ltrain,1)-size(feat_out,1)+1,size(Ltrain,1)];
        end
    end
end

% step 2) feature classification
labels=Ltrain;
feats=Ftrain;
labels_wsi=Lwsi;

KFolds=5;
ACC=[];           % save the final acc
CC=zeros(3,3);    % save the confusion matrix
ItNum=5;         % number of iterations for cross validations
ComNum=50;        % number of PCA components
PCA_usage=1;      % whether use PCA


for cc=1:ItNum
    Indices=crossvalind('Kfold',labels_wsi,KFolds);
    validations=labels_wsi;
    
    Allpredictions=zeros(size(feats,1),1);
    for k=1:KFolds
        OriTest=(Indices==k);
        OriTrain=~OriTest;
        
        OriIndex=false(size(feats,1),1);
        ind=find(OriTrain==1);
        for s=1:length(ind)
            sn=ind(s);
            ss=Lse(sn,1);
            se=Lse(sn,2);
            OriIndex(ss:se)=1;
        end
        
        trainingPredictors=feats(OriIndex,:);
        trainingResponse=labels(OriIndex);
        weigs=Ww(OriIndex);
        
        testingPredictors=feats(~OriIndex,:);
        testingResponse=labels(~OriIndex);
        
        
        % predict for each selected tile
        if PCA_usage==1
            %         [pcaCoefficients,pcaScores,~,~,explained,pcaCenters]=pca(trainingPredictors);  %%??
            %         explainedVarianceToKeepAsFraction=95/100;
            %         numComponentsToKeep = find(cumsum(explained)/sum(explained) >= explainedVarianceToKeepAsFraction, 1);
            %         pcaCoefficients = pcaCoefficients(:,1:numComponentsToKeep);
            %         trainingPredictors2=pcaScores(:,1:numComponentsToKeep);                                                                %%??
            
            [pcaCoefficients, pcaScores, ~, ~, explained, pcaCenters] = pca(...
                trainingPredictors,'NumComponents', ComNum);
            trainingPredictors=pcaScores(:,:);
        end
        
        template = templateSVM(...
            'KernelFunction', 'gaussian', ...
            'PolynomialOrder', [], ...
            'KernelScale', 'auto', ...
            'BoxConstraint', 1, ...
            'Standardize', true);
        classificationSVM = fitcecoc(...
            trainingPredictors, ...
            trainingResponse, ...
            'Weights',weigs,...
            'Learners', template, ...
            'Coding', 'onevsone', ...
            'ClassNames', [1; 2; 3]);
        
        if PCA_usage==1
            pcaTransformationFcn=@(x)(x-pcaCenters)*pcaCoefficients;
            svmPredictFcn=@(x)predict(classificationSVM,x);
            validationPredictionFcn=@(x)svmPredictFcn(pcaTransformationFcn(x));
        else
            svmPredictFcn=@(x)predict(classificationSVM,x);
            validationPredictionFcn=@(x)svmPredictFcn(x);
        end
        
        
        [foldPrediction,foldScores]=validationPredictionFcn(testingPredictors);
        
        Allpredictions(~OriIndex)=foldPrediction;
        
        
        
        
        
    end
    % combine for whole slide prediction
    
    for s=1:length(labels_wsi)
        ss=Lse(s,1);
        se=Lse(s,2);
        %validations(s)=mode(Allpredictions(ss:se));
        cnn=zeros(1,3);
        for cc=1:3
            cnn(cc)=sum((Allpredictions(ss:se)==cc).*Ww(ss:se));
        end
        [mv,md]=max(cnn);
        validations(s)=md;
    end
    
    
    correctPredictions = (validations == labels_wsi);
    ACC=[ACC;sum(correctPredictions)/numel(labels_wsi)];
    
    C = confusionmat(labels_wsi,validations); %% confusion matrxi computation
    CC=CC+C;
end
CC=CC/ItNum;
accl=CC(1,1)/sum(CC(1,:))
acci=CC(2,2)/sum(CC(2,:))
acch=CC(3,3)/sum(CC(3,:))
accm=mean(ACC);

% CC=round(CC);
% save(strcat('C:\Users\xuh3\Desktop\prostate-project\journal_features\plots\','proposed_confusion_SVM.mat'),'CC');




