%% -----main function-----
% the main function to process .svs image for saving AP clustering selected
% image patches
% author: Hongming Xu, Ph.D., Cleveland Clinic
% feel free to use it, but should give original author credits
% email: mxu@ualberta.ca
% code tested on Win10 system
% requriment: matlab-openslide toolbox freely available online

close all;clc;
addpath(genpath('Y:\projects\xhm_code_repos\matlab_repository\toolboxes\openslide-matlab-master\'));
addpath(genpath('..\utility_funcs\'))

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
imagePath={'E:\data\blca_mutationBurden\blca_wsi\'}; % 362 well-quality patients
%imagePath={'E:\data\blca_mutationBurden\blca_wsi\',
%           'E:\data\blca_mutationBurden\blca_wsi2\'}; % 24 not very good quality images

% only used for saving debuging outputs, you can freely change it to
% another direcctory
debugOutput='E:\Hongming\projects\tcga-bladder-mutationburden\history_analysis\debug_output\';

% outputs from this function, you can freely change them based on your
% applications
clusterImgOutput={'E:\Hongming\projects\tcga-bladder-mutationburden\tiles_output\P_CN_20X\'};   % save selected tumor tiles
apfreqOutput='E:\Hongming\projects\tcga-bladder-mutationburden\feature_output\P_CN_20X\apfreq\'; % save representative tumor tile name, and the number of
        % patches in that class

%0) load saved svm tumor detection classifier
load(strcat('../step01)_tumor_versus_nontumor/','SVM_cubic_model.mat'));

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
        
        feat_tumor_cluster=[feat_tumor,(top_left_tumor+bottom_right_tumor)/2];
        
        %% history tests
        % %         [idx,indf]=feature_clustering(feat_tumor_cluster,50); % k-means clustering
        % %         freq=hist(idx,1:length(unique(idx)));
        
        [idx,indf]=feature_clustering_AP(feat_tumor_cluster);       % AP clustering
        freq=zeros(1,length(indf));
        for t=1:length(indf)
            freq(t)=sum(idx==indf(t));
        end
        
        if debug==1
            %% for generating figures shown in the paper
            mpdc10 = distinguishable_colors(length(indf));
            temp=top_left(logical(ylabel),:);
            figure,imshow(RGB)
            for i=1:size(temp,1)
                tl=temp(i,:);
                ind=find(idx(i)==indf);
                pos=[tl(2) tl(1) tileSize(2)-1 tileSize(1)-1];
                hold on, rectangle('Position',pos,'EdgeColor','g','LineWidth',2);
                hold on,plot(tl(2)+tileSize(2)/2,tl(1)+tileSize(1)/2,'*','Color',mpdc10(ind,:),'MarkerSize',15, 'LineWidth',10);
            end
            
            temp2=temp(indf,:);
            
            for t=1:size(temp2,1)
                tp=temp2(t,:);
                hold on,plot(tp(2)+tileSize(2)/2,tp(1)+tileSize(1)/2,'r*','MarkerSize',15, 'LineWidth',10);
            end
            
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
        
        top_left_np=top_left_tumor(indf,:);
        bottom_right_np=bottom_right_tumor(indf,:);
        %%top_left_np=top_left_tumor;
        %%bottom_right_np=bottom_right_tumor;
        
        % 8) optionall for computing lbp feautres - history tests
        if LBP_baseline==true
            feat40=feat_tumor(indf,:);
            if any(mag==mag10)
                levelforRead=find(mag==mag10,1);
                feat72=xu_textureComputation_II(top_left_np,bottom_right_np,slidePtr,levelforRead,mag10,magCoarse);
            else
                magToUseAbove=min(mag(mag>mag10));
                levelforRead=find(mag==magToUseAbove);
                feat72=xu_textureComputation_II(top_left_np,bottom_right_np,slidePtr,levelforRead,mag10,magCoarse,magToUseAbove);
            end
            
            feat_out=[feat40,feat72,freq'];
            featureOutput={'E:\Hongming\projects\tcga-bladder-mutationburden\feature_output\P_CN_20X\1)lbp\'};
            %save(strcat(featureOutput{ip},imgs(im).name,'.mat'),'feat_out');
        end
        
        % 9) save tumor tiles based on ap clustering for subsequent feature extraction 
        if any(mag==mag10)
            levelforRead=find(mag==mag10,1);
            feat_cell=xu_save_selected_image_patches(top_left_np,bottom_right_np,slidePtr,levelforRead,mag10,magCoarse,imgs(im).name(1:23),clusterImgOutput{ip},freq',HE);
        else
            magToUseAbove=min(mag(mag>mag10));
            levelforRead=find(mag==magToUseAbove);
            feat_cell=xu_save_selected_image_patches(top_left_np,bottom_right_np,slidePtr,levelforRead,mag10,magCoarse,imgs(im).name(1:23),clusterImgOutput{ip},freq',HE,magToUseAbove);
        end
        
        % 10) save representative tumor tile name, and the number of
        % patches in that class
        save(strcat(apfreqOutput,imgs(im).name,'.mat'),'feat_cell');
        
    end
end


