

%% -----main function-----
% generate tumor maps for overlapping tils+tmb+tumor together

% author: Hongming Xu, Ph.D., Cleveland Clinic
% feel free to use it, but should give original author credits
% email: mxu@ualberta.ca
% code tested on Win10 system
% requriment: matlab-openslide toolbox freely available online

close all;clc;
addpath(genpath('Y:\projects\xhm_code_repos\matlab_repository\toolboxes\openslide-matlab-master\'));
addpath(genpath('../../utility_funcs/'))

% ---- main parameer settings -----%
magCoarse=2.5;
magFine=5;
mag10=20;
ppt=0.8;             % above this threshold is selected as segmented regions
thrWhite=210;        % above this threshold is considered as whilte background pixels
tileSize=[256,256]./2;

debug=0;             % if 1: to observer intermediate pictures
LBP_baseline=false;  % if false not computing LBP baseline texture features for comparison
% ----- end parameter settings ----%

% ----- input-output-path settings ---%

% tcga_blca .svs slide data can also be accessed from hwang lab shared
% space or downloading from tcga data portal
%imagePath={'E:\data\blca_mutationBurden\blca_wsi\'}; % 362 well-quality patients
imagePath={'E:\data\blca_mutationBurden\blca_wsi\',...
           'E:\data\blca_mutationBurden\blca_wsi2\'}; % 24 not very good quality images

% used for saving outputs maps 
Output_path='../heatmap_blca/tumor_maps/';


%0) load saved svm tumor detection classifier
load(strcat('../../step01)_tumor_versus_nontumor/','SVM_cubic_model.mat'));

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
        if debug==1
            xu_debugShownTiles(RGB,bwTissue,top_left,tileSize);
            %saveas(gcf,strcat(debugOutput,imgs(im).name,'.jpg'));
            close all;
        end
        
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
        
        tumor_pred=zeros(size(bwTissue));
        for t=1:size(top_left_tumor,1)
            tumor_pred(top_left_tumor(t,1):bottom_right_tumor(t,1),top_left_tumor(t,2):bottom_right_tumor(t,2))=1;
        end
        tumor_pred2=imresize(tumor_pred,0.125/2,'nearest');
        imwrite(tumor_pred2*255,strcat(Output_path,imgs(im).name(1:23),'.png'))
        
        
        %% testing recommented by the reviewer -- clustering using features compuated at 10x
        % %         if any(mag==10)
        % %             levelforRead=find(mag==10,1);
        % %             feat_tumor=xu_textureComputation(top_left_tumor,bottom_right_tumor,slidePtr,levelforRead,10,magCoarse);
        % %         else
        % %             magToUseAbove=min(mag(mag>10));
        % %             levelforRead=find(mag==magToUseAbove);
        % %             feat_tumor=xu_textureComputation(top_left_tumor,bottom_right_tumor,slidePtr,levelforRead,10,magCoarse,magToUseAbove);
        % %         end
        %% end--testing
        
    end
end


