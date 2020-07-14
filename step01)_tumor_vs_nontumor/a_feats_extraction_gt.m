%% Author: Hongming Xu, Postdoc, Hwang Lab 2018
% function purpose: 
% input:
% -- xmlfolder path: annotation files
% -- imagePath: pathology slides
% output: save .mat file features into the folder: featureOutput
% 
% 1) read mannualy labeled regions (e.g., tumor) from .svs file,
% 2) generate small patches, compute and save features for training a tumor classifer
% V0.0: test for TCGA_BLCA pathology slides
% 3) email me if you cannot run the program: mxu@ualberta.ca

% NOTE THAT: you can freely change our program to process your data slides
%            but if you want to process our data, you must contact authours
%            for data access permission, we cannot upload data online due
%            ot its size

clear vars;

% you need to download matlab-openlside package
addpath(genpath('Y:\projects\xhm_code_repos\matlab_repository\toolboxes\openslide-matlab-master\'));            % add openslide package
%addpath(genpath('Y:\projects\xhm_code_repos\matlab_repository\my_codes\'));                                     % add reqiured functions
addpath(genpath('..\utility_funcs\'))

%% Z disk: hwangt-share network disk (for our own dataset)
xmlfolder={'Z:\Datasets\Pathology_Slides\TCGA_BLCA_WSI_TumorGT\tumorWithLabel\',...                      % image path with annotated .xml files
           'Z:\Datasets\Pathology_Slides\TCGA_BLCA_WSI_TumorGT\nontumorWithLabel\'};
imagePath={'Z:\Datasets\Pathology_Slides\TCGA_BLCA_WSI_TumorGT\tumorWithLabel\',...                      % image path with .svs format slides
           'Z:\Datasets\Pathology_Slides\TCGA_BLCA_WSI_TumorGT\nontumorWithLabel\'};

% -- parameter settings -- %
magCoarse=2.5;
magFine=5.0;
pp=0.5;             % above this threshold is selected as segmented regions
ppt=0.5;            % tissue foureground should be above this threshold
tileSize=[128,128]; % corresponds to magCoarse resolution
thrWhite=230;

debug=0;
%debugOutput={'E:\Hongming\projects\tcga-bladder-mutationburden\debugoutput\tumor\',...
%    'E:\Hongming\projects\tcga-bladder-mutationburden\debugoutput\nontumor\'};

featureOutput={'.\tumor_prediction_features\tumor\',...
               '.\tumor_prediction_features\nontumor\'};

% -- end parameter settings --%

count=zeros(length(imagePath),1); % count the number of image tiles
temp=[];                          % count the number of image tiles
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
            PatientID=imgs(im).name(1:23);
            Patname=strcat(PatientID,'.mat');
            save(strcat(featureOutput{ip},Patname),'feat');
        end
    end
end


