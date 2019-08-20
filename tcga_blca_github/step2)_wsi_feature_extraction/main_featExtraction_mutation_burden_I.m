%% main function
%% this is the version used for bioinformatics paper
% the main function to process .svs image for save AP clustering selected
% image patches
% author: Hongming Xu, Ph.D., Cleveland Clinic
% feel free to use it, but should give original author credits

close all;clc;
addpath(genpath('E:\matlab_repository\toolboxes\openslide-matlab-master\'));
addpath(genpath('E:\matlab_repository\misc\'));

% ------------- for image color-normaliztion, due to cmputational
% complexity not used
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
% ------------ end for image color normalization

magCoarse=2.5;
magFine=5;
mag10=20;
mag20=20;
ppt=0.8;         %above this threshold is selected as segmented regions
debug=1;

thrWhite=210;
numc=50;
tileSize=[256,256]./2;

%imagePath={'E:\blca_mutationBurden\blca_wsi\'}; % 362 well-quality
%patients
imagePath={'E:\blca_mutationBurden\tumor_detection\'}; % 24 not very good quality images

debugOutput='E:\Hongming\projects\tcga-bladder-mutationburden\debug_output\';

featureOutput={'E:\Hongming\projects\tcga-bladder-mutationburden\feature_output\10)P_E_CN\1)lbp\'};
clusterImgOutput={'E:\Hongming\projects\tcga-bladder-mutationburden\tiles_output\P_E_CN\'};
apfreqOutput='E:\Hongming\projects\tcga-bladder-mutationburden\feature_output\10)P_E_CN\apfreq\';

%0) load saved svm tumor detection classifier
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
            %close all;
        end
        
        %5) feature computation
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
        
        % 6) feature clustering
        feat_tumor=feat(logical(ylabel),:);
        top_left_tumor=top_left(logical(ylabel),:);
        bottom_right_tumor=bottom_right(logical(ylabel),:);
        
        feat_tumor_cluster=[feat_tumor,(top_left_tumor+bottom_right_tumor)/2];
        
        
        %[idx,indf]=feature_clustering(feat_tumor_cluster,numc); % k-means clustering
        %freq=hist(idx,1:length(unique(idx)));
        
        
        [idx,indf]=feature_clustering_AP(feat_tumor_cluster);    % AP clustering
        freq=zeros(1,length(indf));
        for t=1:length(indf)
            freq(t)=sum(idx==indf(t));
        end
        
        
        
        feat40=feat_tumor(indf,:);
        
        if debug==1
            % for generating figures
%             mpdc10 = distinguishable_colors(length(indf));
%             temp=top_left(logical(ylabel),:);
%             figure,imshow(RGB)
%             for i=1:size(temp,1)
%                 tl=temp(i,:);
%                 ind=find(idx(i)==indf);
%                 pos=[tl(2) tl(1) tileSize(2)-1 tileSize(1)-1];
%                 hold on, rectangle('Position',pos,'EdgeColor','g','LineWidth',2);
%                 hold on,plot(tl(2)+tileSize(2)/2,tl(1)+tileSize(1)/2,'*','Color',mpdc10(ind,:),'MarkerSize',15, 'LineWidth',10);
%             end
            
%             temp2=temp(indf,:);
%             
%             for t=1:size(temp2,1)
%                 tp=temp2(t,:);
%                 hold on,plot(tp(2)+tileSize(2)/2,tp(1)+tileSize(1)/2,'r*','MarkerSize',15, 'LineWidth',10);
%             end
            
            
            temp=top_left(logical(ylabel),:);
            temp2=temp(indf,:);
            for t=1:size(temp2,1)
                tp=temp2(t,:);
                hold on,plot(tp(2)+tileSize(2)/2,tp(1)+tileSize(1)/2,'r*','MarkerSize',15);
                hold on,text(tp(2)+tileSize(2)/2,tp(1)+tileSize(1)/2,num2str(idx(t)));
            end
            %saveas(gcf,strcat(debugOutput,imgs(im).name,'.jpg'));
            close all;
        end
        
        
        
        % 7) extract features at higher magnification
        top_left_np=top_left_tumor(indf,:);
        bottom_right_np=bottom_right_tumor(indf,:);
        %%top_left_np=top_left_tumor;
        %%bottom_right_np=bottom_right_tumor;
        if any(mag==mag10)
            levelforRead=find(mag==mag10,1);
            feat72=xu_textureComputation_II(top_left_np,bottom_right_np,slidePtr,levelforRead,mag10,magCoarse);
        else
            magToUseAbove=min(mag(mag>mag10));
            levelforRead=find(mag==magToUseAbove);
            feat72=xu_textureComputation_II(top_left_np,bottom_right_np,slidePtr,levelforRead,mag10,magCoarse,magToUseAbove);
        end
        
        feat_out=[feat40,feat72,freq'];
        save(strcat(featureOutput{ip},imgs(im).name,'.mat'),'feat_out');
        
        %saveas(gcf,strcat(debugOutput{ip},imgs(im).name,'.jpg'));
        %close all;
        
        
        % 8) extract cellular features
        if any(mag==mag10)
            levelforRead=find(mag==mag10,1);
            
            feat_cell=xu_cellular_feature_test(top_left_np,bottom_right_np,slidePtr,levelforRead,mag10,magCoarse,imgs(im).name(1:23),clusterImgOutput{ip},freq',HE);
            
        else
            magToUseAbove=min(mag(mag>mag10));
            levelforRead=find(mag==magToUseAbove);
            
            feat_cell=xu_cellular_feature_test(top_left_np,bottom_right_np,slidePtr,levelforRead,mag10,magCoarse,imgs(im).name(1:23),clusterImgOutput{ip},freq',HE,magToUseAbove);
            
        end
        
        save(strcat(apfreqOutput,imgs(im).name,'.mat'),'feat_cell');
        %feat_out=[feat72];
        %save(strcat(featureOutput{ip},imgs(im).name,'.mat'),'feat_out');
        
        
    end
end


