% this function is used for 5-folds cross valition by SVM classifier with linear kernel
% we test it for 3 class classifications g6, g7 and >=g8
% author: Hongming Xu, Cleveland Clinic, May 2018
% you may modify or use it, but you need give the credit for original
% author


function SVM_classifier_simple_mean_AllComparison

clear vars;

featPath={'E:\Hongming\projects\tcga-bladder-mutationburden\feature_output\10)P_E_APC\2)xception\'};


Ftrain=[];
Ltrain=[];

load('E:\Hongming\projects\tcga-bladder-mutationburden\Hongming_codes\blca_MutBurdens.mat');  % TCGA TMB values
%load('E:\Hongming\projects\tcga-bladder-mutationburden\Hongming_codes\blca_MutBurdens_Yunku.mat'); % Yunku provided values
%load('E:\Hongming\projects\tcga-bladder-mutationburden\Hongming_codes\blca_MutBurdensII.mat'); % TCGA use the second division criterion
% step 1) feature organization
patID=[];
ii=1;

for ip=1:length(featPath)
    fp=featPath{ip};
    
    mats=dir(fullfile(fp,'*.mat'));
    
    for im=1:numel(mats)
        matfilename=strcat(fp,mats(im).name);
        load(matfilename);
        
        patientID=mats(im).name(1:23);

       
               
        
        feat_mean=mean(image_features);
        
        temp=strcmp(blca_MutBurdens(:,1),patientID(1:12));
        
        catg=blca_MutBurdens(temp==1,3);
        if strcmp(catg,'Low')
            Llow=1*ones(size(feat_mean,1),1);
            Ltrain=[Ltrain;Llow];
            Ftrain=[Ftrain;feat_mean];
            patID{ii}=patientID;
            ii=ii+1;
        elseif strcmp(catg,'Mid')
            Lmid=3*ones(size(feat_mean,1),1);
            %Ltrain=[Ltrain;Lmid];
            %Ftrain=[Ftrain;feat_mean];
        elseif strcmp(catg,'High')
            Lhigh=3*ones(size(feat_mean,1),1);
            Ltrain=[Ltrain;Lhigh];
            Ftrain=[Ftrain;feat_mean];
            patID{ii}=patientID;
            ii=ii+1;
        else
            disp('wrong!!!!!!');
        end
        
    end
end

FT=table(Ftrain,Ltrain);
FT.Properties.VariableNames={'features','classes'};
%
%saveName=strcat(featmatoutput,'FT.mat');
%save(saveName,'FT')


% step 2) feature classification
labels=Ltrain;
feats=Ftrain;

KFolds=10;
ACC=[];           % save the final acc
CC=zeros(2,2);    % save the confusion matrix
ItNum=1;         % number of iterations for cross validations
ComNum=100;        % number of PCA components
PCA_usage=1;      % whether use PCA

validm=[];
SSC=zeros(length(labels),2);
pred50=[];
SPE=[];
SEN=[];
for cc=1:ItNum
    
    %Indices=crossvalind('Kfold',numel(labels),KFolds);
    Indices=randperm(numel(labels));
    
    validations=labels;
    validationScores=NaN(length(labels),2);
    %for k=1:KFolds
    for k=1:numel(labels)
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
            'ClassNames', [1; 3]);
        
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
        validationScores(OriTest,:)=foldScores;
    end
    
    validm=[validm,validations];
    correctPredictions = (validations == labels);
    ACC=[ACC;sum(correctPredictions)/numel(labels)];
    
    
    specificity=sum(labels(labels==1)==validations(labels==1))/(sum(labels==1));
    sensitivity=sum(labels(labels==3)==validations(labels==3))/(sum(labels==3));
    
    SPE=[SPE;specificity];
    SEN=[SEN;sensitivity];
    
    C = confusionmat(labels,validations); %% confusion matrxi computation
    CC=CC+C;
    
    SSC=SSC+validationScores;
    pred50=[pred50,validations];
end
CC=CC/ItNum;
SSC=SSC./ItNum;
accl=CC(1,1)/sum(CC(1,:))
acci=CC(2,2)/sum(CC(2,:))
%acch=CC(3,3)/sum(CC(3,:))
accm=mean(ACC);
spem=mean(SPE);
senm=mean(SEN);

pred=ones(length(labels),1);
pred(SSC(:,2)>SSC(:,1))=3;

pred_bladder=cell(length(pred),3);
pred_bladder(:,1)=patID;
pred_bladder(pred==1,2)={'Low'};
pred_bladder(pred==3,2)={'High'};
pred_bladder(labels==1,3)={'Low'}
pred_bladder(labels==3,3)={'High'}

%save(strcat('E:\Hongming\projects\tcga-bladder-mutationburden\Hongming_codes\step7)_plot_figures\','score_label2.mat'),'labels','SSC');
% CC=round(CC);
% save(strcat('C:\Users\xuh3\Desktop\prostate-project\journal_features\plots\','proposed_confusion_SVM.mat'),'CC');




