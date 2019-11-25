%%----read me----%%
% in this version, i) we try to select 10 or 20 patches from detected tumor
% regions
% ii) the patches selection is based on nucleid densities (i.e., from highest)



%% main function for selecting patches for deep learning training

close all;clc;
addpath(genpath('E:\matlab_repository\toolboxes\openslide-matlab-master\'));
addpath(genpath('E:\matlab_repository\misc\'));

magCoarse=2.5;
%magMid=5.0;
magFine=5;
mag10=10;
mag20=20;
ppt=0.8; % above this threshold is selected as segmented regions
%thrOverStain=0.25;
debug=0;

thrWhite=210;
numc=10;
tileSize=[256,256]./2;
np=10;


imagePath={'E:\blca_mutationBurden\low\','E:\blca_mutationBurden\mid\','E:\blca_mutationBurden\high\'};
debugOutput={'E:\Hongming\projects\tcga-bladder-mutationburden\debugoutput\3)low\','E:\Hongming\projects\tcga-bladder-mutationburden\debugoutput\3)mid\','E:\Hongming\projects\tcga-bladder-mutationburden\debugoutput\3)high\'};

featureOutput={'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\4)mutation_predictionII\low\','E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\4)mutation_predictionII\mid\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\4)mutation_predictionII\high\'};

clusterImgOutput={'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_imagesII\low\','E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_imagesII\mid\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_imagesII\high\'};

load(strcat('E:\Hongming\projects\tcga-bladder-mutationburden\Hongming_codes\step1)_tumor_versus_nontumor\','SVM_cubic_model.mat'));

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
        
        %2) pre-processing to get binary masks
        [bwTissue]=wsi_preprocess_tissue(RGB,thrWhite,tileSize(1)*tileSize(2));
        
        %3) get image tile locations
        [top_left,bottom_right]=xu_SelectImageTiles_VII(bwTissue,ppt,tileSize);
        
        
        if debug==0
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
        ylabel=trainedModel_SVM_cubic.predictFcn(ff);
        
        
        temp=top_left(logical(ylabel),:);
        for t=1:size(temp,1)
            tp=temp(t,:);
            hold on,plot(tp(2)+tileSize(2)/2,tp(1)+tileSize(1)/2,'r*','MarkerSize',15);
        end
        
        
        % 6) patche selection
        top_left_tumor=top_left(logical(ylabel),:);
        bottom_right_tumor=bottom_right(logical(ylabel),:);
        feat_tumor=feat(logical(ylabel),:);               % tumor patch feature at 5x
        
        
        if any(mag==mag10)
            levelforRead=find(mag==mag10,1);
            [top_left_np,bottom_right_np,ind]=xu_patch_selection(np,top_left_tumor,bottom_right_tumor,slidePtr,levelforRead,mag10,magCoarse);
        else
            magToUseAbove=min(mag(mag>mag10));
            levelforRead=find(mag==magToUseAbove);
            [top_left_np,bottom_right_np,ind]=xu_patch_selection(np,top_left_tumor,bottom_right_tumor,slidePtr,levelforRead,mag10,magCoarse,magToUseAbove);
        end
        feat40=feat_tumor(ind,:);
        
        for t=1:size(top_left_np,1)
            tp=top_left_np(t,:);
            hold on,plot(tp(2)+tileSize(2)/2,tp(1)+tileSize(1)/2,'y*','MarkerSize',12);
        end
        saveas(gcf,strcat(debugOutput{ip},imgs(im).name,'.jpg'));
        close all;
        
        % 7) extract features at higher magnification
        if any(mag==mag10)
            levelforRead=find(mag==mag10,1);
            feat_np=xu_textureComputation_II(top_left_np,bottom_right_np,slidePtr,levelforRead,mag10,magCoarse);
        else
            magToUseAbove=min(mag(mag>mag10));
            levelforRead=find(mag==magToUseAbove);
            feat_np=xu_textureComputation_II(top_left_np,bottom_right_np,slidePtr,levelforRead,mag10,magCoarse,magToUseAbove);
        end
        
        
        feat_out=[feat40,feat_np];
        save(strcat(featureOutput{ip},imgs(im).name,'.mat'),'feat_out');
        
        
        
        % 8) extract cellular features
        if any(mag==mag20)
            levelforRead=find(mag==mag20,1);
            feat_cell=xu_cellular_feature(top_left_np,bottom_right_np,slidePtr,levelforRead,mag20,magCoarse,imgs(im).name(1:23),clusterImgOutput{ip});
        else
            magToUseAbove=min(mag(mag>mag20));
            levelforRead=find(mag==magToUseAbove);
            feat_cell=xu_cellular_feature(top_left_np,bottom_right_np,slidePtr,levelforRead,mag20,magCoarse,imgs(im).name(1:23),clusterImgOutput{ip},magToUseAbove);
        end
        
        
        
    end
end


