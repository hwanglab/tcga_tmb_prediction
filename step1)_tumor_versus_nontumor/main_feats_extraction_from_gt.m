%% Author: Hongming Xu, Postdoc, Hwang Lab 2018
% functions: 
% 1) read mannualy labeled regions (e.g., tumor) from .svs file,
% 2) generate small patches, compute and save features for training a tumor classifer
% V0.0: test and use for TCGA bladder pathology slides

clear vars;

addpath(genpath('E:\matlab_repository\toolboxes\openslide-matlab-master\'));            % add openslide package
addpath(genpath('E:\matlab_repository\misc\'));                                         % add reqiured functions
addpath(genpath('E:\Hongming\projects\tcga-bladder-mutationburden\Hongming_codes\'));   % add reqiured functions

xmlfolder={'Z:\Datasets\TCGA_BLCA_WSI_TumorGT\tumorWithLabel\',...                      % image path with annotated .xml files
    'Z:\Datasets\TCGA_BLCA_WSI_TumorGT\nontumorWithLabel\'};
imagePath={'Z:\Datasets\TCGA_BLCA_WSI_TumorGT\tumorWithLabel\',...                      % image path with .svs format slides
    'Z:\Datasets\TCGA_BLCA_WSI_TumorGT\nontumorWithLabel\'};

magCoarse=2.5;
magFine=5.0;
pp=0.5;             % above this threshold is selected as segmented regions
ppt=0.5;            % tissue foureground should be above this threshold
tileSize=[128,128]; % corresponds to magCoarse resolution
thrWhite=230;

debug=1;
debugOutput={'E:\Hongming\projects\tcga-bladder-mutationburden\debugoutput\tumor\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\debugoutput\nontumor\'};
featureOutput={'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\1)tumor_prediction\tumor\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\1)tumor_prediction\nontumor\'};

count=zeros(length(imagePath),1);
temp=[];
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
        
        
        if exist(fullFileName,'file')
            %1) read magCoarse image
            RGB=wsi_read(slidePtr,objectivePower,downsampleFactors,width,height,magCoarse);
            
            %2) read GT binary mask
            [m,n]=size(RGB(:,:,1));
            bwGT=wsi_read_GT(fullFileName,m,n,tt);
            %figure,imshow(bwGT);
            
            %3) select GT tiles
            
            if ip==1  % select tumor regions
                gimg=rgb2gray(RGB);
                bwTissue=(gimg<=thrWhite);
            else      % select non-tumor regions
                bwTissue=ones(m,n);
            end
            [top_left,bottom_right]=xu_SelectImageTiles_VI(bwTissue,bwGT,pp,ppt,tileSize);
            
            count(ip,1)=size(top_left,1)+count(ip,1);
            temp=[temp,size(top_left,1)];
            fprintf('number of patches=%d\n',size(top_left,1));
            
            if debug==1
                xu_debugShownTiles(RGB,bwGT,top_left,tileSize);
                temp01=bwperim(bwGT);
                [r,c]=find(temp01);
                hold on,plot(c,r,'r*'); hold off;
                %saveas(gcf,strcat(debugOutput{ip},imgs(im).name,'.jpg'));
                % save debugged figs
                close all;
            end
            
            %4) feature computation and storation
            if any(mag==magFine)
                levelforRead=find(mag==magFine,1);
                feat=xu_textureComputation(top_left,bottom_right,slidePtr,levelforRead,magFine,magCoarse);
            else
                magToUseAbove=min(mag(mag>magFine));
                levelforRead=find(mag==magToUseAbove);
                feat=xu_textureComputation(top_left,bottom_right,slidePtr,levelforRead,magFine,magCoarse,magToUseAbove);
            end
             

            % save computed features into .mat file for training tumor
            % detection classifier
%             PatientID=imgs(im).name(1:23);
%             Patname=strcat(PatientID,'.mat');
%             save(strcat(featureOutput{ip},Patname),'feat');
        end
    end
end
t=0;

