%%-- A Simple Example for BLCA TMB Prediction
% step1: select and save representative tumor tiles

% --- input: whole slide image .svs
% --- output1: saved selected tumor tiles
% --- output2: a .mat file saving selected tiles' name and the number of
% tiles within that cluster (AP Clustering)

% author: Hongming Xu, Ph.D., Cleveland Clinic
% feel free to use it, but should give original author credits
% email: mxu@ualberta.ca
% code tested on Win10 system
% requriment: matlab-openslide toolbox freely available online

close all;clc;
addpath(genpath('.\openslide-matlab-master\'));
addpath(genpath('..\utility_funcs\'))

% ---- main parameer settings -----%
magCoarse=2.5;       % low magnification to process the whole tissue
magFine=5;           % magnification for lbp feature exraction
mag10=20;            % high magnifiction for saving selected tiles
ppt=0.8;             % above this threshold is selected as segmented regions
thrWhite=210;        % above this threshold is considered as whilte background pixels
tileSize=[256,256]./2; % tile size at magCoarse magnification
output='./step1_output/';
% ----- end parameter settings ----%

%0) load saved svm tumor detection classifier
load(strcat('../step01)_tumor_versus_nontumor/','SVM_cubic_model.mat'));  % for simplicity, used svm tumor detecor here

%% start wsi processing
img_name='TCGA-BL-A3JM-01Z-00-DX1.33E53972-CEA4-4D84-A5D2-7DAD7B0C27F8.svs'; % you need download this sample wsi form tcga and put it in this folder

file1=fullfile('./',img_name);
slidePtr=openslide_open(file1);
[mppX,mppY,width,height,numberOfLevels,...
    downsampleFactors,objectivePower]=openslide_get_slide_properties(slidePtr);
mag=objectivePower./round(downsampleFactors);

%1) read magCoarse image
RGB=wsi_read(slidePtr,objectivePower,downsampleFactors,width,height,magCoarse);

%2) pre-processing to get binary masks with tissue regions as the
%foreground
[bwTissue]=wsi_preprocess_tissue(RGB,thrWhite,tileSize(1)*tileSize(2));

%3) obtain the H,E stain vectors used for color normalization later
HE=extract_hestains(RGB,bwTissue);

%4) get image tile locations
[top_left,bottom_right]=xu_SelectImageTiles_VII(bwTissue,ppt,tileSize);

% for showing intermediate results

xu_debugShownTiles(RGB,bwTissue,top_left,tileSize);
%saveas(gcf,strcat(debugOutput,imgs(im).name,'.jpg'));
close all;


%5) feature computation at low magnification (5x) for ap clustering
if any(mag==magFine)
    levelforRead=find(mag==magFine,1);
    feat=xu_textureComputation(top_left,bottom_right,slidePtr,levelforRead,magFine,magCoarse);
else
    magToUseAbove=min(mag(mag>magFine));
    levelforRead=find(mag==magToUseAbove);
    feat=xu_textureComputation(top_left,bottom_right,slidePtr,levelforRead,magFine,magCoarse,magToUseAbove);
end

% 6) tumor patch classification
ff=table(feat);
ff.Properties.VariableNames={'features'};
[ylabel,scores]=trainedModel_SVM_cubic.predictFcn(ff);

% 7) feature clustering (AP clustering)
feat_tumor=feat(logical(ylabel),:);   %% tumore tile lbp features

top_left_tumor=top_left(logical(ylabel),:);
bottom_right_tumor=bottom_right(logical(ylabel),:);

feat_tumor_cluster=[feat_tumor,(top_left_tumor+bottom_right_tumor)/2];

% 8) ap clustering
[idx,indf]=feature_clustering_AP(feat_tumor_cluster);       % AP clustering
freq=zeros(1,length(indf));
for t=1:length(indf)
    freq(t)=sum(idx==indf(t));
end


top_left_np=top_left_tumor(indf,:);
bottom_right_np=bottom_right_tumor(indf,:);

% 9) save tumor tiles based on ap clustering for subsequent feature extraction
if any(mag==mag10)
    levelforRead=find(mag==mag10,1);
    feat_cell=xu_save_selected_image_patches(top_left_np,bottom_right_np,slidePtr,levelforRead,mag10,magCoarse,img_name(1:23),output,freq',HE);
else
    magToUseAbove=min(mag(mag>mag10));
    levelforRead=find(mag==magToUseAbove);
    feat_cell=xu_save_selected_image_patches(top_left_np,bottom_right_np,slidePtr,levelforRead,mag10,magCoarse,img_name(1:23),output,freq',HE,magToUseAbove);
end

% 10) save representative tumor tile name, and the number of
% patches in that class
save(strcat(output,img_name(1:23),'.mat'),'feat_cell');




