% purpose:
%- 1) generate tmb predicton heatmaps for visualization of intermediate tmb patients
%- 2) save tmb prediction probabilities for entropy computation (main purpose)

% By Hongming Xu, CCF, 2019
% email: mxu@ualberta.ca
% reqirements: two tool boxex are reqiures, you can freely download them
% online
% -1: matlab-openslide
% -2: cbrewer

close all;clc;
addpath(genpath('Y:\projects\xhm_code_repos\matlab_repository\toolboxes\openslide-matlab-master\'));
addpath(genpath('Y:\projects\xhm_code_repos\matlab_repository\toolboxes\cbrewer\'));
addpath(addpath(genpath('..\utility_funcs\')))

magCoarse=2.5;
magFine=5;
ppt=0.8;            % above this threshold is selected as segmented regions
thrWhite=210;       % above this threshold considered as background
tileSize=[256,256]./2;

debug=0;
savePredMap=1;
showHeatmap=1;

%--- path settings ---%
% tcga_blca .svs slides
imagePath={'E:\data\blca_mutationBurden\blca_wsi\'}; % 362 well-quality patients
%imagePath={'E:\blca_mutationBurden\blca_wsi2\'}; % 24 not very good quality images

debugOutput='E:\Hongming\projects\tcga-bladder-mutationburden\debug_output\';

% the path to save prediction scores for tumor tiles
score_path='.\prediction_scores_mid\';
% path to output heatmaps for visulizationa and entropy computations
heatmap_output='.\heatmap_blca\color_maps\mid\';
tmb_output='.\heatmap_blca\mat_files\mid\';
%--- end path settings ---%

load(strcat('..\step01)_tumor_versus_nontumor\','SVM_cubic_model.mat'));
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
        
        mats_cc=dir(fullfile(score_path,strcat(imgs(im).name(1:23),'*.mat')));
        
        if length(mats_cc)==1
            load(strcat(score_path,mats_cc.name));
            
            heat_cc=dir(fullfile(tmb_output,strcat(imgs(im).name(1:23),'*.mat')));
            %if isempty(heat_cc) % if not got it in the first time
            
            %1) read magCoarse image
            RGB=wsi_read(slidePtr,objectivePower,downsampleFactors,width,height,magCoarse);
            
            %2) pre-processing to get binary masks
            [bwTissue]=wsi_preprocess_tissue(RGB,thrWhite,tileSize(1)*tileSize(2));
            
            %3) get image tile locations
            [top_left,bottom_right]=xu_SelectImageTiles_VII(bwTissue,ppt,tileSize);
            
            if debug==1
                xu_debugShownTiles(RGB,bwTissue,top_left,tileSize);
                close all;
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
            
            % 6) feature clustering
            feat_tumor=feat(logical(ylabel),:);
            top_left_tumor=top_left(logical(ylabel),:);
            bottom_right_tumor=bottom_right(logical(ylabel),:);
            
            feat_tumor_cluster=[feat_tumor,(top_left_tumor+bottom_right_tumor)/2];
            
            [idx,indf]=feature_clustering_AP(feat_tumor_cluster);    % AP clustering
            freq=zeros(1,length(indf));
            for t=1:length(indf)
                freq(t)=sum(idx==indf(t));
            end
            
            if length(indf)==size(foldScores,1)
                pp_map=NaN(size(bwTissue));
                for i=1:length(idx)
                    cc=idx(i);
                    ind=find(cc==indf);
                    %if ind<=length(foldScores)
                    ss=foldScores(ind,2);
                    tl=top_left_tumor(i,:);
                    br=bottom_right_tumor(i,:);
                    pp_map(tl(1):br(1),tl(2):br(2))=ss;
                    %end
                end
                if savePredMap==1
                    tmb_map=pp_map;
                    save(strcat(tmb_output,imgs(im).name(1:23),'.mat'),'tmb_map','top_left_tumor','bottom_right_tumor');
                end
                
                if showHeatmap==1
                    % use our found matlab functions
                    pp_map2=imresize(pp_map,0.125,'nearest');
                    %pp_map2(pp_map2<0.01)=NaN;
                    CT=cbrewer('div', 'RdYlBu', 64);
                    colormap(flipud(CT));
                    heatmap_v0(pp_map2,[],[],[], 'NaNColor', [1 1 1]);
                    
                    cbh = colorbar('southoutside') ; %Create Colorbar
                    cbh.Ticks = linspace(min(min(pp_map2)), max(max(pp_map2)), 2) ; %Create 8 ticks from zero to 1
                    temp=cell(1,2);
                    temp{1}='Low TMB';
                    temp{2}='High TMB';
                    cbh.TickLabels = temp;
                    saveas(gcf,strcat(heatmap_output,imgs(im).name,'.jpg'));
                    close all;
                end
            else
                fprintf('filename=%s\n%d--no matching!!!!!!!!!!!!',file1,im);
            end
        end
        %end
    end
end


