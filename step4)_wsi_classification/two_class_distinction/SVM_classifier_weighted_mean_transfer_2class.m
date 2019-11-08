% this function is used for 5-folds cross valition by SVM classifier with linear kernel
% we test it for 3 class classifications g6, g7 and >=g8
% author: Hongming Xu, Cleveland Clinic, May 2018
% you may modify or use it, but you need give the credit for original
% author


function SVM_classifier_weighted_mean_transfer_2class

clear vars;
%featPath={'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\2)mutation_prediction\low\',...
%          'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\2)mutation_prediction\mid\',...
%          'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\2)mutation_prediction\high\'};

featPath_cluster={'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\7)mutation_prediction_clustering\low\',...
                  'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\7)mutation_prediction_clustering\mid\',...
                  'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\7)mutation_prediction_clustering\high\'};

% featPath={'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\7)transfer_learning_k50\vgg19\low\',...
%           'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\7)transfer_learning_k50\vgg19\mid\',...
%           'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\7)transfer_learning_k50\vgg19\high\'};

featPath={'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\7)transfer_learning_k50\inceptionv3\low\',...
          'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\7)transfer_learning_k50\inceptionv3\mid\',...
          'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\7)transfer_learning_k50\inceptionv3\high\'}; %inceptionV3

%featmatoutput='E:\Hongming\projects\tcga-bladder-mutationburden\Hongming_codes\step4)_wsi_classification\';

Ftrain=[];
Ltrain=[];

% step 1) feature organization
pat=cell(363,1);
ii=1;

for ip=1:length(featPath)
    fp=featPath{ip};
    cp=featPath_cluster{ip};
    mats=dir(fullfile(fp,'*.mat'));
    
    for im=1:numel(mats)
        matfilename=strcat(fp,mats(im).name);
        load(matfilename);
        
        patientID=mats(im).name(1:23);
        mats_cc=dir(fullfile(cp,strcat(mats(im).name(1:23),'*.mat')));
        if length(mats_cc)==1
            load(strcat(cp,mats_cc.name));
            coeff=feat_out(:,end)./sum(feat_out(:,end));
            coeff00=repmat(coeff,1,size(image_features,2));
            feat_mean=sum(image_features.*coeff00,1);
            
%             feat00=feat_out(:,1:end-1);
%             coeff00=repmat(coeff,1,size(feat00,2));
%             feat_mean01=sum(feat00.*coeff00,1);
            
            feat_mean=[feat_mean];
            
            %Ftrain=[Ftrain;feat_mean];
            
            if ip==1
                Llow=ones(size(feat_mean,1),1);
                Ltrain=[Ltrain;Llow];
                Ftrain=[Ftrain;feat_mean];
            elseif ip==2
                Lmid=2*ones(size(feat_mean,1),1);
                %Ltrain=[Ltrain;Lmid];
                %Ftrain=[Ftrain;feat_mean];
            else
                Lhigh=2*ones(size(feat_mean,1),1);
                Ltrain=[Ltrain;Lhigh];
                Ftrain=[Ftrain;feat_mean];
            end
        else
            disp('impossible~~~~~~');
            break;
        end
        
        
    end
end

% FT=table(Ftrain,Ltrain);
% FT.Properties.VariableNames={'features','classes'};
% 
% saveName=strcat(featmatoutput,'FT.mat');
% save(saveName,'FT')


% step 2) feature classification
labels=Ltrain;
feats=Ftrain;

KFolds=5;
ACC=[];           % save the final acc
CC=zeros(2,2);    % save the confusion matrix
ItNum=50;         % number of iterations for cross validations
ComNum=30;        % number of PCA components
PCA_usage=1;      % whether use PCA

validm=[];
ks='auto';
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
            'KernelScale', ks, ...
            'BoxConstraint', 1, ...
            'Standardize', true);
        
%         template = templateSVM(...
%         'KernelFunction', 'polynomial', ...
%         'PolynomialOrder', 3, ...
%         'KernelScale', ks, ...
%         'BoxConstraint', 1, ...
%         'Standardize', true);
    
        classificationSVM = fitcecoc(...
            trainingPredictors, ...
            trainingResponse, ...
            'Learners', template, ...
            'Coding', 'onevsone', ...
            'ClassNames', [1; 2]);
        
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
accl=CC(1,1)/sum(CC(1,:));
acci=CC(2,2)/sum(CC(2,:));
%acch=CC(3,3)/sum(CC(3,:))
accm=mean(ACC);

% CC=round(CC);
% save(strcat('C:\Users\xuh3\Desktop\prostate-project\journal_features\plots\','proposed_confusion_SVM.mat'),'CC');



