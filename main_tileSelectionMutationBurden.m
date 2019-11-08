%% main function this is for dr.Lee's work to extract LBP features

close all;clc;
addpath(genpath('E:\matlab_repository\toolboxes\openslide-matlab-master\'));
addpath(genpath('E:\matlab_repository\misc\'));

magCoarse=2.5;
magMid=10;
magFine=20;
pp=0.01;             %% choose blocks with nuclei pixels above than 10% of top value ???? make sure patches have nuclei
ppt=0.8;             %% choose blocks with tissue pixels above than 90% of block size ???? make sure pathes mainly have tissues

debug=0;

thrWhite=220;
thrOverStain=0.4*255;
ss=10;
Tumor_tileSize=[256,256]./2;
Tumor_augmentStep=0; %% augmenting steps to achieve half overlapping
mnum=20;
mapping=getmapping(16,'riu2');
mapping2=getmapping(8,'riu2');
testing=0;

imagePath={'E:\bladder_project_DX\dataset\'};


featureOutput={'E:\bladder_project\featOutput_tcga\'};

debugOutput='E:\bladder_project\debug_output\';

for ip=1:length(imagePath)
    if ~testing
        imgPath=imagePath{ip};
    else
        imgPath=testimagePath{ip};
    end
    imgs=dir(fullfile(imgPath,'*.svs'));
    
    for im=10:numel(imgs)
        file1=fullfile(imgPath,imgs(im).name);
        fprintf('filename=%s\n%d',file1,im);
        slidePtr=openslide_open(file1);
        [mppX,mppY,width,height,numberOfLevels,...
            downsampleFactors,objectivePower]=openslide_get_slide_properties(slidePtr);
        mag=objectivePower./round(downsampleFactors);
        
        
        %1) read magCoarse image
        RGB=wsi_read(slidePtr,objectivePower,downsampleFactors,width,height,magCoarse);
        
        %2) pre-processing to get binary masks
        [bwTissue,bwHimg]=wsi_preprocess_bladder_forLee(RGB,thrWhite,thrOverStain,round(Tumor_tileSize(1)*Tumor_tileSize(2)*5));
        
        %3) get image tile locations
        [top_left,bottom_right]=xu_SelectImageTiles_V(bwTissue,bwHimg,pp,ppt,Tumor_tileSize,mnum);
        %[top_left,bottom_right]=xu_SelectImageTiles_IV(bwTissue,bwHimg,pp,ppt,Tumor_tileSize,Tumor_augmentStep,mnum); % step 4: select interested image tiles this is for 2.5x by default
        
       
        
        
        
        if debug==1
            xu_debugShownTiles(RGB,bwTissue,top_left,Tumor_tileSize);
            saveas(gcf,strcat(debugOutput,imgs(im).name,'.jpg'));
            close all;
        end
        t=0;
        %4) feature computation
        if any(mag==magMid)
            levelforRead=find(mag==magMid,1);
            PatF=xu_texturalFeats(top_left,bottom_right,slidePtr,levelforRead,magMid,magCoarse,mapping,mapping2);
        else
            magToUseAbove=min(mag(mag>magMid));
            levelforRead=find(mag==magToUseAbove);
            PatF=xu_texturalFeats(top_left,bottom_right,slidePtr,levelforRead,magMid,magCoarse,mapping,mapping2,magToUseAbove);
        end
        
        temp=strsplit(featureOutput{ip},'\');
        %Patname=strcat(char(temp(end-1)),'-',num2str(im),'.mat');
        PatientID=imgs(im).name(1:16);
        Patname=strcat(PatientID,'.mat');
        save(strcat(featureOutput{ip},Patname),'PatF');

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
%             saveas(gcf,strcat(debugOutput,imgs(im).name,'.jpg'));
%             close all;
%         end
%        close all;
    end
    
    t=0;
end


