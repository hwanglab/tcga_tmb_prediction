% this function is used bladder cancer tumor detection evaluation
% reported in the paper
% author: Hongming Xu, Cleveland Clinic, May 2018
% you may modify or use it, but you need give the credit for original
% author


function svm_tumor_detection_paper_results

clear vars;
featureoutput={'E:\Hongming\projects\tcga-bladder-mutationburden\feature_output\1)tumor_prediction\tumor\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\feature_output\1)tumor_prediction\nontumor\'};


% step 1) feature organization
Ftrain=[];
Ltrain=[];

temp2=[];
for ip=1:length(featureoutput)
    featpath=featureoutput{ip};
    mats=dir(fullfile(featpath,'*.mat'));
    
    for im=1:numel(mats)
        matfilename=strcat(featpath,mats(im).name);
        %fprintf('load filename=%s\n%d',mats(im).name,im);
        load(matfilename);
        Ftrain=[Ftrain;feat];
        temp2=[temp2,size(feat,1)];
        
        if ip==1 % tumor
            Ltumor=ones(size(feat,1),1);
            Ltrain=[Ltrain;Ltumor];
        else
            Lnontumor=zeros(size(feat,1),1);
            Ltrain=[Ltrain;Lnontumor];
        end
    end
end

FT=table(Ftrain,Ltrain);
FT.Properties.VariableNames={'features','classes'};


% step 2) feature classification
labels=Ltrain;
feats=Ftrain;

KFolds=5;
ItNum=1;          % number of iterations for cross validations
ComNum=20;        % number of PCA components
PCA_usage=0;      % whether use PCA

for cc=1:ItNum
    
    Indices=crossvalind('Kfold',numel(labels),KFolds);
    %validations=labels;
    %validationScores=NaN(length(labels),2);
    for k=1:1    % 20% of patches for validation evaluation
   
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
            %         trainingPredictors2=pcaScores(:,1:numComponentsToKeep);                                             
            
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
            'ClassNames', [0; 1]);
        
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
        
        
        %validations(OriTest)=foldPrediction;
        %validationScores(OriTest,:)=foldScores;
    end
    
    % prediction accuracy, specificity and sensitivity reproted in the
    % paper
    ACC=sum(testingResponse==foldPrediction)/numel(testingResponse);
    specificity=sum(testingResponse(testingResponse==0)==foldPrediction(testingResponse==0))/(sum(testingResponse==0));
    sensitivity=sum(testingResponse(testingResponse==1)==foldPrediction(testingResponse==1))/(sum(testingResponse==1))
    
    %validm=[validm,validations];
    %correctPredictions = (validations == labels);
    %ACC=[ACC;sum(correctPredictions)/numel(labels)];
    
end

% save prediction scores for ploting ROC curves
%save(strcat('E:\Hongming\projects\tcga-bladder-mutationburden\Hongming_codes\step1)_tumor_versus_nontumor\','bladder_score_label.mat'),'foldScores','testingResponse');





