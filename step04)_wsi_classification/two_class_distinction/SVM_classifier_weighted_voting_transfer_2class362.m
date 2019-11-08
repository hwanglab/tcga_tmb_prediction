% this function is used for 5-folds cross valition by SVM classifier with linear kernel
% we test it for 3 class classifications g6, g7 and >=g8
% author: Hongming Xu, Cleveland Clinic, May 2018
% you may modify or use it, but you need give the credit for original
% author


function SVM_classifier_weighted_voting_transfer_2class362

clear vars;

ap_cluster={'E:\Hongming\projects\tcga-bladder-mutationburden\feature_output\10)norm_test_20x\apfreq\'};
featPath={'E:\Hongming\projects\tcga-bladder-mutationburden\feature_output\10)norm_test_20x\5)inceptionv3\'};

Ftrain=[];
Ltrain=[];
Lwsi=[];
Lse=[];
Ww=[];

load('E:\Hongming\projects\tcga-bladder-mutationburden\Hongming_codes\blca_MutBurdens.mat');  % TCGA TMB values

% step 1) feature organization
for ip=1:length(featPath)
    fp=featPath{ip};
    cp=ap_cluster{ip};
    mats=dir(fullfile(fp,'*.mat'));
    
    for im=1:numel(mats)
        matfilename=strcat(fp,mats(im).name);
        load(matfilename);
        
        patientID=mats(im).name(1:23);
        mats_cc=dir(fullfile(cp,strcat(mats(im).name(1:23),'*.mat')));
        if length(mats_cc)==1
            load(strcat(cp,mats_cc.name));
            
            apn=cell2mat(feat_cell(:,2));
            apn=apn./sum(apn);
            %coeff=apn./sum(apn);
            
            
            temp=strcmp(blca_MutBurdens(:,1),patientID(1:12));
            
            catg=blca_MutBurdens(temp==1,3);
            if strcmp(catg,'Low')
                Llow=ones(size(image_features,1),1);
                Ltrain=[Ltrain;Llow];
                Lwsi=[Lwsi;1];
                Lse=[Lse;size(Ltrain,1)-size(image_features,1)+1,size(Ltrain,1)];
                
                Ww=[Ww;apn];
                Ftrain=[Ftrain;image_features];
                
            elseif strcmp(catg,'Mid')
%                 Lmid=2*ones(size(image_features,1),1);
%                 Ltrain=[Ltrain;Lmid];
%                 Lwsi=[Lwsi;2];
%                 Lse=[Lse;size(Ltrain,1)-size(image_features,1)+1,size(Ltrain,1)];
%                 
%                 Ww=[Ww;apn];
%                 Ftrain=[Ftrain;image_features];
            elseif strcmp(catg,'High') 
                Lhigh=2*ones(size(image_features,1),1);
                Ltrain=[Ltrain;Lhigh];
                Lwsi=[Lwsi;2];
                Lse=[Lse;size(Ltrain,1)-size(image_features,1)+1,size(Ltrain,1)];
                
                Ww=[Ww;apn];
                Ftrain=[Ftrain;image_features];
            end
        else
            disp('impossible~~~~~~');
            break;
        end
    end
end

% step 2) feature classification
labels=Ltrain;
feats=Ftrain;
labels_wsi=Lwsi;

KFolds=5;
ACC=[];           % save the final acc
CC=zeros(2,2);    % save the confusion matrix
ItNum=10;         % number of iterations for cross validations
ComNum=100;        % number of PCA components
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
            'ClassNames', [1; 2]);
        
%         classificationSVM = fitcecoc(...
%             trainingPredictors, ...
%             trainingResponse, ...
%             'Learners', template, ...
%             'Coding', 'onevsone', ...
%             'ClassNames', [1; 2]);
        
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
        cnn=zeros(1,2);
        for cc=1:2
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
%acch=CC(3,3)/sum(CC(3,:))
accm=mean(ACC);

% CC=round(CC);
% save(strcat('C:\Users\xuh3\Desktop\prostate-project\journal_features\plots\','proposed_confusion_SVM.mat'),'CC');




