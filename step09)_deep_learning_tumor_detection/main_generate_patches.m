%% Author: Hongming Xu, Postdoc, Hwang Lab 2018
% function usage: read mannualy labeled regions (e.g., tumor) from .svs
% file, and generate small patches, compute and save features for training a tumor classifer
% V0.0: test and use for TCGA bladder pathology slides

clear vars;

addpath(genpath('E:\matlab_repository\toolboxes\openslide-matlab-master\'));
addpath(genpath('E:\matlab_repository\misc\'));
addpath(genpath('E:\Hongming\projects\tcga-bladder-mutationburden\Hongming_codes\'));

xmlfolder={'Z:\Datasets\TCGA_BLCA_WSI_TumorGT\refinement_version\tumorWithLabel\'};
imagePath={'Z:\Datasets\TCGA_BLCA_WSI_TumorGT\refinement_version\tumorWithLabel\'};

magCoarse=2.5;
magFine=5.0;
mag10=20;
mag20=20;
pp=0.5; % above this threshold is selected as segmented regions
ppt=0.2; % tissue foureground should above this threshold
ppn=0.5;
tileSize=[128,128]./2; % corresponds to magCoarse resolution
thrWhite=200;
wid=150;

debug=0;
debugOutput={'E:\Hongming\projects\tcga-bladder-mutationburden\debugoutput\tumor\','E:\Hongming\projects\tcga-bladder-mutationburden\debugoutput\nontumor\'};
%patchOutput={'Y:\DL_tumor_detection\data_bladder\train\tumor\','Y:\DL_tumor_detection\data_bladder\train\non_tumor\'};  % first all output into the train folders
patchOutput={'Y:\DL_tumor_detection\trail01\data_bladder_20x\train\tumor\','Y:\DL_tumor_detection\trail01\data_bladder_20x\train\non_tumor\'};  % first all output into the train folders
countTumor=0;
countNonTumor=0;
tempT=[];
tempN=[];
for ip=1:length(imagePath)
    
    imgPath=imagePath{ip};
    xmlPath=xmlfolder{ip};
    
    imgs=dir(fullfile(imgPath,'*.svs'));
    
    
    for im=1:numel(imgs)
        file1=fullfile(imgPath,imgs(im).name);
        fprintf('filename=%s\n%d',file1,im);
        slidePtr=openslide_open(file1);
        [mppX,mppY,width,height,numberOfLevels,...
            downsampleFactors,objectivePower]=openslide_get_slide_properties(slidePtr);
        mag=objectivePower./round(downsampleFactors);
        
        
        tt=objectivePower/magCoarse;
        xmlname=strcat(imgs(im).name(1:end-3),'xml');
        fullFileName=fullfile(xmlPath,xmlname);
        
        %1) read magCoarse image
        RGB=wsi_read(slidePtr,objectivePower,downsampleFactors,width,height,magCoarse);
        
        %         figure,imshow(RGB);
        %         gimg=rgb2gray(RGB);
        %         bwTissue=(gimg<=thrWhite);
        %         figure,imshow(bwTissue);
        % [r,c]=find(bwTissue);
        % hold on,plot(c,r,'b.');
        %       close all;
        
        if exist(fullFileName,'file')
            
            %2) read GT binary mask
            [m,n]=size(RGB(:,:,1));
            bwGT=wsi_read_GT(fullFileName,m,n,tt);
            
            bwGT_nonTumor0=~bwGT;
            bwGT_nonTumor=false(size(bwGT_nonTumor0));
            bwGT_nonTumor(wid:size(bwGT_nonTumor0,1)-wid,wid:size(bwGT_nonTumor0,2)-wid)=bwGT_nonTumor0(wid:size(bwGT_nonTumor0,1)-wid,wid:size(bwGT_nonTumor0,2)-wid);
            %figure,imshow(bwGT);
            
            %3) select GT tiles
            gimg=rgb2gray(RGB);
            bwTissue=(gimg<=thrWhite);
            [top_left,bottom_right]=xu_SelectImageTiles_VIII(bwTissue,bwGT,pp,ppt,tileSize);
            
            [top_left2,bottom_right2]=xu_SelectImageTiles_VIII(bwTissue,bwGT_nonTumor,pp,ppn,tileSize);
            
            countTumor=size(top_left,1)+countTumor;
            countNonTumor=size(top_left2,1)+countNonTumor;
            
            tempT=[tempT,size(top_left,1)];
            tempN=[tempN,size(top_left2,1)];
            
            
            if debug==1
                xu_debugShownTiles(RGB,bwGT,top_left,tileSize);
                temp01=bwperim(bwGT);
                [r,c]=find(temp01);
                hold on,plot(c,r,'r*'); hold off;
                saveas(gcf,strcat(debugOutput{2},imgs(im).name,'.jpg'));
                close all;
            end
            
            
           
            if any(mag==mag10)
                levelforRead=find(mag==mag10,1);
                xu_generate_maps(top_left,bottom_right,slidePtr,levelforRead,mag10,magCoarse,imgs(im).name(1:23),patchOutput{1});
                
                xu_generate_maps(top_left2,bottom_right2,slidePtr,levelforRead,mag10,magCoarse,imgs(im).name(1:23),patchOutput{2});
            else
                magToUseAbove=min(mag(mag>mag10));
                levelforRead=find(mag==magToUseAbove);
                xu_generate_maps(top_left,bottom_right,slidePtr,levelforRead,mag10,magCoarse,imgs(im).name(1:23),patchOutput{1},magToUseAbove);
                
                xu_generate_maps(top_left2,bottom_right2,slidePtr,levelforRead,mag10,magCoarse,imgs(im).name(1:23),patchOutput{2},magToUseAbove)
                
            end
            
        else
            bwGT_nonTumor0=true(size(RGB,1),size(RGB,2));
            bwGT_nonTumor=false(size(bwGT_nonTumor0));
            bwGT_nonTumor(wid:size(bwGT_nonTumor0,1)-wid,wid:size(bwGT_nonTumor0,2)-wid)=bwGT_nonTumor0(wid:size(bwGT_nonTumor0,1)-wid,wid:size(bwGT_nonTumor0,2)-wid);
            
            gimg=rgb2gray(RGB);
            bwTissue=(gimg<=thrWhite);
            [top_left2,bottom_right2]=xu_SelectImageTiles_VIII(bwTissue,bwGT_nonTumor,pp,ppn,tileSize);
            
            
            countNonTumor=size(top_left2,1)+countNonTumor;
            tempN=[tempN,size(top_left2,1)];
            
%             if debug==1
%                 xu_debugShownTiles(RGB,bwGT_nonTumor,top_left2,tileSize);
%                 temp01=bwperim(bwGT_nonTumor);
%                 [r,c]=find(temp01);
%                 hold on,plot(c,r,'r*'); hold off;
%                 saveas(gcf,strcat(debugOutput{2},imgs(im).name,'.jpg'));
%                 close all;
%             end
            
            if any(mag==mag10)
                levelforRead=find(mag==mag10,1);
                
                xu_generate_maps(top_left2,bottom_right2,slidePtr,levelforRead,mag10,magCoarse,imgs(im).name(1:23),patchOutput{2});
            else
                magToUseAbove=min(mag(mag>mag10));
                levelforRead=find(mag==magToUseAbove);

                xu_generate_maps(top_left2,bottom_right2,slidePtr,levelforRead,mag10,magCoarse,imgs(im).name(1:23),patchOutput{2},magToUseAbove)
                
            end
            
        end
        
    end
end
t=0;

