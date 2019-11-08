%% main function
%% this is the version used for bioinformatic paper

close all;clc;
addpath(genpath('E:\matlab_repository\toolboxes\openslide-matlab-master\'));
addpath(genpath('E:\matlab_repository\misc\'));
addpath(genpath('E:\matlab_repository\toolboxes\cbrewer\'));

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


magCoarse=2.5;
magFine=5;
mag10=20;
mag20=20;
ppt=0.8; % above this threshold is selected as segmented regions
debug=0;

thrWhite=210;
numc=50;

tileSize=[256,256]./2;

imagePath={'E:\blca_mutationBurden\blca_wsi\'}; % 362 well-quality
%patients
%imagePath={'E:\blca_mutationBurden\tumor_detection\'}; % 24 not very good quality images

debugOutput='E:\Hongming\projects\tcga-bladder-mutationburden\debug_output\';

featureOutput={'E:\Hongming\projects\tcga-bladder-mutationburden\feature_output\10)P_E_CN\1)lbp\'};
clusterImgOutput={'E:\Hongming\projects\tcga-bladder-mutationburden\tiles_output\P_E_CN\'};
apfreqOutput='E:\Hongming\projects\tcga-bladder-mutationburden\feature_output\10)P_E_CN\apfreq\';

load(strcat('E:\Hongming\projects\tcga-bladder-mutationburden\Hongming_codes\step1)_tumor_versus_nontumor\','SVM_cubic_model.mat'));

score_path='E:\Hongming\projects\tcga-bladder-mutationburden\Hongming_codes\step4)_wsi_classification\two_class_distinction\prediction_scores\';
%heatmap_output='E:\Hongming\projects\tcga-bladder-mutationburden\heatmap_output\';
tmb_output='E:\Hongming\projects\tcga-bladder-mutationburden\tmb_blca\';


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
        
        if length(mats_cc)==1  % if bladder patient belonging to 253 selected TMB prediction groups
            load(strcat(score_path,mats_cc.name));
            
            heat_cc=dir(fullfile(tmb_output,strcat(imgs(im).name(1:23),'*.mat')));
            if isempty(heat_cc) % if not got it in the first time
                
                %1) read magCoarse image
                RGB=wsi_read(slidePtr,objectivePower,downsampleFactors,width,height,magCoarse);
                
                %2) pre-processing to get binary masks
                [bwTissue]=wsi_preprocess_tissue(RGB,thrWhite,tileSize(1)*tileSize(2));
                
                HE=extract_hestains(RGB,bwTissue);
                
                %3) get image tile locations
                [top_left,bottom_right]=xu_SelectImageTiles_VII(bwTissue,ppt,tileSize);
                
                
                if debug==1
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
                    tmb_map=pp_map;
                    save(strcat(tmb_output,imgs(im).name(1:23),'.mat'),'tmb_map','top_left_tumor','bottom_right_tumor');
                    
%                     % use our found matlab functions
%                     pp_map2=imresize(pp_map,0.125,'nearest');
%                     %pp_map2(pp_map2<0.01)=NaN;
%                     CT=cbrewer('div', 'RdYlBu', 64);
%                     colormap(flipud(CT));
%                     heatmap_v0(pp_map2,[],[],[], 'NaNColor', [1 1 1]);
%                     
%                     cbh = colorbar('southoutside') ; %Create Colorbar
%                     cbh.Ticks = linspace(min(min(pp_map2)), max(max(pp_map2)), 2) ; %Create 8 ticks from zero to 1
%                     temp=cell(1,2);
%                     temp{1}='Low TMB';
%                     temp{2}='High TMB';
%                     cbh.TickLabels = temp;
%                     saveas(gcf,strcat(heatmap_output,imgs(im).name,'_uncertain.jpg'));
%                     close all;
                else
                    fprintf('filename=%s\n%d--no matching!!!!!!!!!!!!',file1,im);
                end
            end
            % use toolbox plotly
            %pp_map2=flipud(pp_map2);
            %data={struct('z',pp_map2,'showscale','false','type','heatmap')};
            %response = plotly(data, struct('filename', 'basic-heatmap', 'fileopt', 'overwrite'));
            %plot_url = response.url;
            %heatmap(pp_map2,[],[],[],'ColorMap', @cool, 'NaNColor', [1 1 1], 'colorbar', true);
            
            
            %             feat40=feat_tumor(indf,:);
            %
            %             if debug==1
            %                 % for generating figures
            %                 %             mpdc10 = distinguishable_colors(length(indf));
            %                 %             temp=top_left(logical(ylabel),:);
            %                 %             figure,imshow(RGB)
            %                 %             for i=1:size(temp,1)
            %                 %                 tl=temp(i,:);
            %                 %                 ind=find(idx(i)==indf);
            %                 %                 pos=[tl(2) tl(1) tileSize(2)-1 tileSize(1)-1];
            %                 %                 hold on, rectangle('Position',pos,'EdgeColor','g','LineWidth',2);
            %                 %                 hold on,plot(tl(2)+tileSize(2)/2,tl(1)+tileSize(1)/2,'*','Color',mpdc10(ind,:),'MarkerSize',15, 'LineWidth',10);
            %                 %             end
            %
            %                 %             temp2=temp(indf,:);
            %                 %
            %                 %             for t=1:size(temp2,1)
            %                 %                 tp=temp2(t,:);
            %                 %                 hold on,plot(tp(2)+tileSize(2)/2,tp(1)+tileSize(1)/2,'r*','MarkerSize',15, 'LineWidth',10);
            %                 %             end
            %
            %
            %                 temp=top_left(logical(ylabel),:);
            %                 temp2=temp(indf,:);
            %                 for t=1:size(temp2,1)
            %                     tp=temp2(t,:);
            %                     hold on,plot(tp(2)+tileSize(2)/2,tp(1)+tileSize(1)/2,'r*','MarkerSize',15);
            %                     hold on,text(tp(2)+tileSize(2)/2,tp(1)+tileSize(1)/2,num2str(idx(t)));
            %                 end
            %saveas(gcf,strcat(debugOutput,imgs(im).name,'.jpg'));
            %close all;
            %end
        end
        
        
        
        
        
    end
end


