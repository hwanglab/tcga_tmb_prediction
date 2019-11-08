% this function is used for 5-folds cross valition by SVM classifier with linear kernel
% we test it for 3 class classifications g6, g7 and >=g8
% author: Hongming Xu, Cleveland Clinic, May 2018
% you may modify or use it, but you need give the credit for original
% author


function cross_validation_SVM_classifier_simple_mean

clear vars;
%featPath={'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\2)mutation_prediction\low\','E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\2)mutation_prediction\mid\',...
%    'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\2)mutation_prediction\high\'};

%featPath={'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\6)mutation_prediction\low\','E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\6)mutation_prediction\mid\',...
%     'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\6)mutation_prediction\high\'};
% 

featPath={'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\7)mutation_prediction_clustering\low\','E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\7)mutation_prediction_clustering\mid\',...
     'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\7)mutation_prediction_clustering\high\'};

Ftrain=[];
Ltrain=[];

% v-2: feature organization for averaging methods
for ip=1:length(featPath)
    fp=featPath{ip};
    mats=dir(fullfile(fp,'*.mat'));
    
    for im=1:numel(mats)
        matfilename=strcat(fp,mats(im).name);
        load(matfilename)
        
        feat_mean=mean(feat_out(:,41:112));
        
        Ftrain=[Ftrain;feat_mean];
        
        if ip==1
            Llow=ones(size(feat_mean,1),1);
            Ltrain=[Ltrain;Llow];
        elseif ip==2
            Lmid=2*ones(size(feat_mean,1),1);
            Ltrain=[Ltrain;Lmid];
        else
            Lhigh=3*ones(size(feat_mean,1),1);
            Ltrain=[Ltrain;Lhigh];
        end
    end
end


% step 2) feature classification
labels=Ltrain;
feats=Ftrain;

KFolds=5;
ACC=[];           % save the final acc
CC=zeros(3,3);    % save the confusion matrix
ItNum=50;         % number of iterations for cross validations
ComNum=30;        % number of PCA components
PCA_usage=1;      % whether use PCA

validm=[];
for cc=1:ItNum
    Indices=crossvalind('Kfold',labels,KFolds);
    validations=labels;
    for k=1:KFolds
        OriTest=(Indices==k);
        OriTrain=~OriTest;
        trainingPredictors=feats(OriTrain,:);
        trainingResponse=labels(OriTrain,:);
        
        
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
        
        
        testingPredictors=feats(OriTest,:);
        testingResponse=labels(OriTest);
        [foldPrediction,foldScores]=validationPredictionFcn(testingPredictors);
        
        
        validations(OriTest)=foldPrediction;
    end
    
    validm=[validm,validations];
    correctPredictions = (validations == labels);
    ACC=[ACC;sum(correctPredictions)/numel(labels)];
    
    C = confusionmat(labels,validations); %% confusion matrxi computation
    CC=CC+C;
end
CC=CC/ItNum;
accl=CC(1,1)/sum(CC(1,:))
acci=CC(2,2)/sum(CC(2,:))
acch=CC(3,3)/sum(CC(3,:))
accm=mean(ACC);

% CC=round(CC);
% save(strcat('C:\Users\xuh3\Desktop\prostate-project\journal_features\plots\','proposed_confusion_SVM.mat'),'CC');




