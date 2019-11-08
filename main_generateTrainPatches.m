%% main function for selecting patches for deep learning training

close all;clc;
addpath(genpath('E:\matlab_repository\toolboxes\openslide-matlab-master\'));
addpath(genpath('E:\matlab_repository\misc\'));

magCoarse=2.5;
magMid=5.0;
magFine=20;
ppt=0.9;             %% choose blocks with tissue pixels above than 90% of block size ???? make sure pathes mainly have tissues

debug=0;

thrWhite=210;


tileSize=[256,256]./2;
mapping=getmapping(16,'riu2');
testing=0;

imagePath={'Z:\Datasets\TCGA_BLCA_WSI_TumorGT\tumorWithLabel\'};
debugOutput='E:\Hongming\projects\tcga-bladder-mutationburden\debugoutput\';

featureOutput={'C:\Users\xuh3\Desktop\projects\tcga-bladder-subtype\featureOutput\luminal_papillary\','C:\Users\xuh3\Desktop\projects\tcga-bladder-subtype\featureOutput\luminal_infiltrated\',...
    'C:\Users\xuh3\Desktop\projects\tcga-bladder-subtype\featureOutput\luminal\','C:\Users\xuh3\Desktop\projects\tcga-bladder-subtype\featureOutput\basal_squamous\','C:\Users\xuh3\Desktop\projects\tcga-bladder-subtype\featureOutput\neuronal\'};

%xmlfolder='Z:\Datasets\TCGA_BLCA_WSI_TumorGT\tumorWithLabel\';   %% tumor path

xmlfolder='Z:\Datasets\TCGA_BLCA_WSI_TumorGT\nontumorWithLabel\';  %% non-tumor path

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
        [bwTissue]=wsi_preprocess_foreground(RGB,thrWhite);
        
        %3) pre-processing to get tumor GT or non-tumor GT masks
        bwGT=wsi_preprocess_GT(xmlfolder,imgs(im).name,objectivePower,magCoarse,size(RGB(:,:,1)));
        
        %3) get image tile locations
        [top_left,bottom_right]=xu_SelectImageTiles_VI(bwTissue,bwGT,ppt,tileSize);
  
        
        if debug==0
            xu_debugShownTiles(RGB,bwTissue,top_left,tileSize);
            saveas(gcf,strcat(debugOutput,imgs(im).name,'.jpg'));
            close all;
        end
        
%         %4) feature computation
%         if any(mag==magMid)
%             levelforRead=find(mag==magMid,1);
%             PatF=xu_texturalFeats(top_left,bottom_right,slidePtr,levelforRead,magMid,magCoarse,mapping,mapping2);
%         else
%             magToUseAbove=min(mag(mag>magMid));
%             levelforRead=find(mag==magToUseAbove);
%             PatF=xu_texturalFeats(top_left,bottom_right,slidePtr,levelforRead,magMid,magCoarse,mapping,mapping2,magToUseAbove);
%         end
%         
%         temp=strsplit(featureOutput{ip},'\');
%         Patname=strcat(char(temp(end-1)),'-',num2str(im),'.mat');
%         PatientID=imgs(im).name(1:16);
%         save(strcat(featureOutput{ip},Patname),'PatF','PatientID');
        
        %        imagename=strcat(num2str(im),'_',imgs(im).name(1:12));
        %         %4) select high resolution image tile
        %         if any(mag==magFine)
        %             levelforRead=find(mag==magFine,1);
        %
        %             xu_generateTiles_II(patchOutput{ip},imagename,bwHimg,top_left,bottom_right,slidePtr,levelforRead,magFine,magCoarse);
        %
        % %             [tf,br]=xu_selectTiles(bwHimg,top_left,bottom_right,slidePtr,levelforRead,magFine,magCoarse);
        % %
        % %             if ~testing
        % %                 xu_generateTiles(patchOutput{ip},imagename,neiN,rotN,ps,tf,br,slidePtr,levelforRead,magFine,magCoarse);
        % %             else
        % %                 xu_generateTiles(testpatchOutput{ip},imagename,neiN,1,ps,tf,br,slidePtr,levelforRead,magFine,magCoarse);     %% no rotation
        % %             end
        %
        %         else
        %             magToUseAbove=min(mag(mag>magFine));
        %             levelforRead=find(mag==magToUseAbove);
        %
        %             xu_generateTiles_II(patchOutput{ip},imagename,bwHimg,top_left,bottom_right,slidePtr,levelforRead,magFine,magCoarse,magToUseAbove);
        %
        % %             [tf,br]=xu_selectTiles(bwHimg,top_left,bottom_right,slidePtr,levelforRead,magFine,magCoarse,magToUseAbove);
        % %
        % %             if ~testing
        % %                 xu_generateTiles(patchOutput{ip},imagename,neiN,rotN,ps,tf,br,slidePtr,levelforRead,magFine,magCoarse,magToUseAbove);
        % %             else
        % %                 xu_generateTiles(testpatchOutput{ip},imagename,neiN,1,ps,tf,br,slidePtr,levelforRead,magFine,magCoarse,magToUseAbove);  %% no rotation
        % %             end
        %         end
        %
        %         if debug==1
        %             tl=tf;
        %             b1=br;
        %             pos=[tl(2) tl(1) b1(2)-tl(2)-1 b1(1)-tl(1)-1];
        %             hold on, rectangle('Position',pos,'EdgeColor','b','LineWidth',6);
        %             hold on,plot(tl(2)+tileSize(2)/2,tl(1)+tileSize(1)/2,'mp','MarkerSize',15,'LineWidth',5);
        %             hold off;
        %             saveas(gcf,strcat(debugOutput,imagename,'.jpg'));
        %             close all;
        %         end
        
    end
end


