%%-- A Simple Example for BLCA TMB Prediction
% step3: classification

% --- input: outputs from step1 and step2
% --- output1: patient-level and tile-level predictions

% author: Hongming Xu, Ph.D., Cleveland Clinic
% feel free to use it, but should give original author credits
% email: mxu@ualberta.ca
% code tested on Win10 system
% requriment: matlab-openslide toolbox freely available online
clc;clear all;
img_id='TCGA-BL-A3JM-01Z-00-DX1';


% 1) load features
load(strcat('./step2_output/',img_id,'_feat.mat'));
% 2) load ap clustering information
load(strcat('./step1_output/',img_id,'.mat'));
apn=cell2mat(feat_cell(:,2));

% weighted mean feature vector
coeff=apn./sum(apn);
coeff00=repmat(coeff,1,size(image_features,2));
feat_mean=sum(image_features.*coeff00,1);

% 3) load tmb_svm_classifier
load('./tmb_svm_classifier.mat');

pcaTransformationFcn=@(x)(x-pcaCenters)*pcaCoefficients;
svmPredictFcn=@(x)predict(ScoreSVMModel,x);
validationPredictionFcn=@(x)svmPredictFcn(pcaTransformationFcn(x));

% 4) patient-level prediction
[pat_Prediction,pat_Scores]=validationPredictionFcn(feat_mean);
fprintf('patient sample %s tmb prediction score is %f\n', img_id, pat_Scores(1,2));
if pat_Scores(1,2)>=0.5
    fprintf('patient sample %s is predicted as tmb HIGH is %f\n', img_id);
else
    fprintf('patient sample %s is predicted as tmb LOW\n', img_id);
end
disp('---------------------------')

% 5) tile-level predictions
[tile_Prediction,tile_Scores]=validationPredictionFcn(image_features);
for i=1:size(image_names,1)
    fprintf('tile %s tmb prediction score is %f\n', image_names(i,:), tile_Scores(i,2));
end
disp('--------------- end ------------')