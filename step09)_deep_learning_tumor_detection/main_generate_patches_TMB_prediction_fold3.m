%% main function for selecting patches for deep learning training

close all;clc;
addpath(genpath('E:\matlab_repository\toolboxes\openslide-matlab-master\'));
addpath(genpath('E:\matlab_repository\misc\'));

% addpath(genpath('E:/matlab_repository/toolboxes/spams-matlab-v2.6-2017-02-27/spams-matlab-v2.6/'));
% addpath(genpath('E:/Hongming/resources/CodeRelease_ColorNormalization-master/SNMF stain separation and color normalization/'));
% targetpath='E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\tissue_images_normalized\TCGA-E2-A14V-01Z-00-DX1.png';
% nstains=2;
% lambda=0.1;  % Use smaller values of the lambda (0.01-0.1) for better reconstruction. however, if the normalized image seems not fine, increase or decrease the value accordingly.
% [Wi, Hi, Hiv]=stainsep(target,nstains,lambda);
% para.Wi=Wi;
% para.Hi=Hi;
% para.Hiv=Hiv;
% para.nstains=nstains;
% para.lambda=lambda;


magCoarse=2.5;
magFine=5;
mag10=20;
mag20=20;
ppt=0.8; % above this threshold is selected as segmented regions
debug=0;

thrWhite=210;

tileSize=[256,256]./2;

%imagePath={'E:\blca_mutationBurden\blca_wsi\'};
imagePath={'E:\blca_mutationBurden\tumor_detection\'};
debugOutput='E:\Hongming\projects\tcga-bladder-mutationburden\debug_output\';

load(strcat('E:\Hongming\projects\tcga-bladder-mutationburden\Hongming_codes\step1)_tumor_versus_nontumor\','SVM_cubic_model.mat'));

%% (1) build a cross validation index
load('E:\Hongming\projects\tcga-bladder-mutationburden\Hongming_codes\blca_MutBurdens.mat');  % TCGA TMB values
wsi_label=[];
for ip=1:length(imagePath)
    
    imgPath=imagePath{ip};
    
    imgs=dir(fullfile(imgPath,'*.svs'));
    
    for im=1:numel(imgs)
        patientID=imgs(im).name(1:23);
        temp=strcmp(blca_MutBurdens(:,1),patientID(1:12));
        
        catg=blca_MutBurdens(temp==1,3);
        if strcmp(catg,'Low')
            wsi_label=[wsi_label;10];
        elseif strcmp(catg,'Mid')
            wsi_label=[wsi_label;20];
        elseif strcmp(catg,'High')
            wsi_label=[wsi_label;30];
        else
            disp('impossible!!!!');
        end
    end
end
rng('default');
rng(1);
KFolds=3;
Indices=crossvalind('Kfold',wsi_label,KFolds);

%% (2) extract & save patches corresponding directories
for ip=1:length(imagePath)
    
    imgPath=imagePath{ip};
    imgs=dir(fullfile(imgPath,'*.svs'));
    
    
    for im=1:numel(imgs)
        file1=fullfile(imgPath,imgs(im).name);
        fprintf('filename=%s\n%d',file1,im);
        slidePtr=openslide_open(file1);
        [mppX,mppY,width,height,numberOfLevels,...
            downsampleFactors,objectivePower]=openslide_get_slide_properties(slidePtr);
        mag=objectivePower./round(downsampleFactors);
        
        patientID=imgs(im).name(1:23);
        temp=strcmp(blca_MutBurdens(:,1),patientID(1:12));
        catg=blca_MutBurdens(temp==1,3);
        if Indices(im)==1 % for valiation   3: fold1, 2: fold2, 1: fold3
            switch catg{1}
                case 'Low'
                    clusterImgOutput='Y:\DL_TMB_prediction\data_bladder_20X_norm\fold03\validation\low\';
                case 'Mid'
                    clusterImgOutput='Y:\DL_TMB_prediction\data_bladder_20X_norm\fold03\validation\mid\';
                case 'High'
                    clusterImgOutput='Y:\DL_TMB_prediction\data_bladder_20X_norm\fold03\validation\high\';
                otherwise
                    disp('impossible!!!!');
            end
        else
            switch catg{1}
                case 'Low'
                    clusterImgOutput='Y:\DL_TMB_prediction\data_bladder_20X_norm\fold03\training\low\';
                case 'Mid'
                    clusterImgOutput='Y:\DL_TMB_prediction\data_bladder_20X_norm\fold03\training\mid\';
                case 'High'
                    clusterImgOutput='Y:\DL_TMB_prediction\data_bladder_20X_norm\fold03\training\high\';
            end
        end
        
        %1) read magCoarse image
        RGB=wsi_read(slidePtr,objectivePower,downsampleFactors,width,height,magCoarse);
        
        %2) pre-processing to get binary masks
        [bwTissue]=wsi_preprocess_tissue(RGB,thrWhite,tileSize(1)*tileSize(2));
        
        HE=extract_hestains(RGB,bwTissue);
        
        %3) get image tile locations
        [top_left,bottom_right]=xu_SelectImageTiles_VII(bwTissue,ppt,tileSize);
        
        
        if debug==1
            xu_debugShownTiles(RGB,bwTissue,top_left,tileSize);
            %saveas(gcf,strcat(debugOutput,imgs(im).name,'.jpg'));
            %close all;
        end
        
        %4) feature computation
        if any(mag==magFine)
            levelforRead=find(mag==magFine,1);
            feat=xu_textureComputation(top_left,bottom_right,slidePtr,levelforRead,magFine,magCoarse);
        else
            magToUseAbove=min(mag(mag>magFine));
            levelforRead=find(mag==magToUseAbove);
            feat=xu_textureComputation(top_left,bottom_right,slidePtr,levelforRead,magFine,magCoarse,magToUseAbove);
        end
        
        % 5) tumor region classification
        ff=table(feat);
        ff.Properties.VariableNames={'features'};
        [ylabel,scores]=trainedModel_SVM_cubic.predictFcn(ff);
        %ylabel=scores(:,2)>0.5;   % use this thresholds to select tumor patches
        
     
        top_left_tumor=top_left(logical(ylabel),:);
        bottom_right_tumor=bottom_right(logical(ylabel),:);
        
           
        if debug==1
            
            temp=top_left(logical(ylabel),:);
            for t=1:size(temp,1)
                tp=temp(t,:);
                hold on,plot(tp(2)+tileSize(2)/2,tp(1)+tileSize(1)/2,'r*','MarkerSize',15);
            end
           
        end
 
        
        
        % 8) extract cellular features
        if any(mag==mag10)
            levelforRead=find(mag==mag10,1);
            
            xu_generate_tiles_forDL(top_left_tumor,bottom_right_tumor,slidePtr,levelforRead,mag10,magCoarse,imgs(im).name(1:23),clusterImgOutput,HE);
            
        else
            magToUseAbove=min(mag(mag>mag10));
            levelforRead=find(mag==magToUseAbove);
            
            xu_generate_tiles_forDL(top_left_tumor,bottom_right_tumor,slidePtr,levelforRead,mag10,magCoarse,imgs(im).name(1:23),clusterImgOutput,HE,magToUseAbove);
            
        end
        
    end
end


