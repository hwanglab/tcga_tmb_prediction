
% main purpose: process mid patients
% - use svm classifier to make predictions on individual selected tumor
% tiles for MID TMB patients
% and save the scores which are used for generating heatmaps
% -- linear svm is trained by using low and high tmb patients
% -- testing is performed on ap selected tumor tiles of mid patients

% author: Hongming Xu, Cleveland Clinic, Jan 2019
% you may modify or use it, but you need give the credit for original
% author
% email: mxu@ualberta.ca


function c_svm_prediction_scores_mid

clear vars;

%% --- these are the indictors which testing we are performing---%
%  --- for more detail experiments see the paper description ---%
techs={'P_E_TD','P_E_CN','P_InceptionV3','P_Resnet50','P_Xception'};
meth=techs{5}; % 1,2,3,4,5 corresponding different testings

% features of tumor tiles output by pretrained deep learning models
switch meth
    case 'P_E_TD'
        ap_cluster={'E:\Hongming\projects\tcga-bladder-mutationburden\feature_output\P_E_TD\apfreq\'};
        featPath={'E:\Hongming\projects\tcga-bladder-mutationburden\feature_output\P_E_TD\3)xception\'};
    case 'P_E_CN'
        ap_cluster={'E:\Hongming\projects\tcga-bladder-mutationburden\feature_output\P_E_CN\apfreq\'};
        featPath={'E:\Hongming\projects\tcga-bladder-mutationburden\feature_output\P_E_CN\2)xception\'};
    case 'P_InceptionV3'
        ap_cluster={'E:\Hongming\projects\tcga-bladder-mutationburden\feature_output\P_CN_20X\apfreq\'};
        featPath={'E:\Hongming\projects\tcga-bladder-mutationburden\feature_output\P_CN_20X\5)inceptionv3\'};
    case 'P_Resnet50'
        ap_cluster={'E:\Hongming\projects\tcga-bladder-mutationburden\feature_output\P_CN_20X\apfreq\'};
        featPath={'E:\Hongming\projects\tcga-bladder-mutationburden\feature_output\P_CN_20X\2)resnet50\'};
    case 'P_Xception'
        ap_cluster={'E:\Hongming\projects\tcga-bladder-mutationburden\feature_output\P_CN_20X\apfreq\'};
        featPath={'E:\Hongming\projects\tcga-bladder-mutationburden\feature_output\P_CN_20X\4)xception\'};
    otherwise
        disp('impossible options!!!!!!!!!!');
end

%% -- feature organization into a table structure---%%
Ftrain=[];
Ltrain=[];

Ftest=[];
Ltest=[];

load('..\blca_MutBurdens.mat');  % TCGA TMB categories by 1-third percentile

patID=[];                        % only save intermediate tmb level patients
ii=1;
tt=[];                           % only used for counting the number of selected tumor tiles

fp=featPath{1};
cp=ap_cluster{1};
mats=dir(fullfile(fp,'*.mat'));

for im=1:numel(mats)
    matfilename=strcat(fp,mats(im).name);
    load(matfilename);
    
    patientID=mats(im).name(1:23);
    mats_cc=dir(fullfile(cp,strcat(mats(im).name(1:23),'*.mat')));
    
    if length(mats_cc)==1
        load(strcat(cp,mats_cc.name));
        
        apn=cell2mat(feat_cell(:,2));
        coeff=apn./sum(apn);
        
        coeff00=repmat(coeff,1,size(image_features,2));
        feat_mean=sum(image_features.*coeff00,1);
        
        temp=strcmp(blca_MutBurdens(:,1),patientID(1:12));
        
        catg=blca_MutBurdens(temp==1,3);
        if strcmp(catg,'Low')
            Llow=1*ones(size(feat_mean,1),1);
            Ltrain=[Ltrain;Llow];
            Ftrain=[Ftrain;feat_mean];
            %patID{ii}=patientID;
            %ii=ii+1;
            tt=[tt,size(image_features,1)];
        elseif strcmp(catg,'Mid')
            Lmid=2*ones(size(feat_mean,1),1);
            Ltest=[Ltest;Lmid];
            Ftest=[Ftest;feat_mean];
            patID{ii}=patientID;
            ii=ii+1;
        elseif strcmp(catg,'High')
            Lhigh=3*ones(size(feat_mean,1),1);
            Ltrain=[Ltrain;Lhigh];
            Ftrain=[Ftrain;feat_mean];
            %patID{ii}=patientID;
            %ii=ii+1;
            tt=[tt,size(image_features,1)];
        else
            disp('wrong!!!!!!');
        end
    else
        disp('impossible~~~~~~');
        break;
    end
end

FT=table(Ftrain,Ltrain);
FT.Properties.VariableNames={'features','classes'};

%saveName=strcat(featmatoutput,'FT.mat');
%save(saveName,'FT')

%% ---SVM training ---%
score_output='.\prediction_scores_mid\';

labels=Ltrain;
feats=Ftrain;

ComNum=100;       % number of PCA components, might be changed for different methods!!!!!!!!!!
PCA_usage=1;      % whether use PCA

rng(100); % for the same indice we could generate

trainingPredictors=feats;
trainingResponse=labels;

% predict for each selected tile
if PCA_usage==1
    
    %% -- PCA by explained variance --%
    %         [pcaCoefficients,pcaScores,~,~,explained,pcaCenters]=pca(trainingPredictors);  %%??
    %         explainedVarianceToKeepAsFraction=95/100;
    %         numComponentsToKeep = find(cumsum(explained)/sum(explained) >= explainedVarianceToKeepAsFraction, 1);
    %         pcaCoefficients = pcaCoefficients(:,1:numComponentsToKeep);
    %         trainingPredictors2=pcaScores(:,1:numComponentsToKeep);                                                                %%??
    
    %% -- PCA by number of components selected --%
    [pcaCoefficients, pcaScores, ~, ~, explained, pcaCenters] = pca(...
        trainingPredictors,'NumComponents', ComNum);
    trainingPredictors=pcaScores(:,:);
end

classificationSVM = fitcsvm(...
    trainingPredictors, ...
    trainingResponse, ...
    'KernelFunction', 'linear', ...
    'PolynomialOrder', [], ...
    'KernelScale', 'auto', ...
    'BoxConstraint', 1, ...
    'Standardize', true, ...
    'ClassNames', [1; 3]);

ScoreSVMModel=fitPosterior(classificationSVM);

if PCA_usage==1
    pcaTransformationFcn=@(x)(x-pcaCenters)*pcaCoefficients;
    svmPredictFcn=@(x)predict(ScoreSVMModel,x);
    validationPredictionFcn=@(x)svmPredictFcn(pcaTransformationFcn(x));
else
    svmPredictFcn=@(x)predict(classificationSVM,x);
    validationPredictionFcn=@(x)svmPredictFcn(x);
end

%2) testing on intermediate patients
save_scores=false;
foldss=[];
pred=[];
for kk=1:length(patID)
    pid=patID{kk};
    matfilename=strcat(fp,pid,'_feat.mat');
    load(matfilename);
    mats_cc=dir(fullfile(cp,strcat(pid,'*.mat')));
    load(strcat(cp,mats_cc.name));
    apn=cell2mat(feat_cell(:,2));
    
    % weighted mean feature vector
    coeff=apn./sum(apn);
    coeff00=repmat(coeff,1,size(image_features,2));
    feat_mean=sum(image_features.*coeff00,1);
    % ending weighted mean feature vector
    
    [foldPrediction,foldScores]=validationPredictionFcn(image_features);
    
    if save_scores==true
        save(strcat(score_output,pid,'.mat'),'foldScores'); % used for
        %generating heatmaps
    end
    
    lowN=sum((foldPrediction==1).*apn);
    highN=sum((foldPrediction==3).*apn);
    if lowN>highN
        foldPrediction2=1;
    else
        foldPrediction2=3;
    end
    
    coeff00=repmat(coeff,1,size(foldScores,2));
    foldScores2=sum(foldScores.*coeff00,1);
    
    pred=[pred;foldPrediction2];
    foldss=[foldss;foldScores2];
end





