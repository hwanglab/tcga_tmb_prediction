function extract_nuclei_feature

imagePath={'E:\blca_mutationBurden\low\','E:\blca_mutationBurden\mid\','E:\blca_mutationBurden\high\'};

seg_path={'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_segs\low\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_segs\int\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_segs\high\'};

image_path={'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_images\low\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_images\mid\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_images\high\'};

debug_path={'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_debugs2\low\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_debugs2\int\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_debugs2\high\'};

feat_output={'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\3)mutation_prediction_cellular_features\low\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\3)mutation_prediction_cellular_features\mid\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\3)mutation_prediction_cellular_features\high\'};


numc=5;  % number of nuclei clusters
method='kmeans';
fb=10;
for ip=1:length(imagePath)
    
    imgPath=imagePath{ip};
    imgs=dir(fullfile(imgPath,'*.svs'));
    
    for im=1:numel(imgs)
        file1=fullfile(imgPath,imgs(im).name);
        fprintf('filename=%s\n%d',file1,im);
        
        patientID=imgs(im).name(1:23);
        imgseg=dir(fullfile(seg_path{ip},strcat(patientID,'*')));
        
        ffn_img=[];
        ffc=cell(numel(imgseg),1);
        for ss=1:numel(imgseg)
            bwn=imread(strcat(seg_path{ip},imgseg(ss).name));
            temp=strsplit(imgseg(ss).name,'.');
            img=imread(strcat(image_path{ip},temp{1},'.png'));
            
            %1) cell nuclei feature extrcation
            fm=regionprops(bwn,'Area','MajorAxisLength','MinorAxisLength','Eccentricity');
            fnc=regionprops(bwn,'Centroid');  % nuclei location
            fnxy=[fnc.Centroid];
            fnxy2=zeros(length(fnc),2);
            fnxy2(:,1)=fnxy(1:2:end);
            fnxy2(:,2)=fnxy(2:2:end);
            ffc{ss,1}=fnxy2;
            
            %             %2) dilated nuclei features
            %             bwd=imdilate(bwn,strel('disk',6));
            %             D=bwdist(~bwd);
            %             D=-D;
            %             D=imimposemin(D,bwn);
            %             Lw=watershed(D);
            %             bwd(Lw==0)=0;
            %             %bwn=xor(bwareaopen(bwn,siz_nuc),bwareaopen(bwn,size_nuc_max));
            %             %bwn=imfill(bwn,'holes');
            
            
            
            %3 nuclei intensity features
            fmm=cell2mat(struct2cell(fm));
            fratio=fmm(2,:)./fmm(3,:);
            
            fc=regionprops(bwn,img(:,:,1),'PixelValues','MeanIntensity');
            fc00=[fc.MeanIntensity];
            fc01=zeros(1,length(fc00));
            for k=1:length(fc00)
                fc01(k)=std(double(fc(k).PixelValues));
            end
            
            ffn=[fmm;fratio;fc00;fc01]'; % six columns are: area, maj, minor, ecc, fratio, mean intensity and entropy
            %ffn=[fmm(1,:);fmm(4,:);fc00]';
            ffn_img=[ffn_img;ffn];
            
            %             img2=imoverlay(img,bwperim(bwn));
            %             figure,imshow(img2);
            %             hold on,plot(fnxy2(:,1),fnxy2(:,2),'g+');
        end
        
%         [idx,C]=feature_clusteringII(ffn_img,numc,method);
%         
%         if strcmp(method,'gaussian')
%             [mmv,L]=max(idx,[],2);
%             L(mmv<0.75)=numc+1;
%         else
%             L=idx;
%         end
%         %         [bv,bi]=sort([sum(L==1),sum(L==2),sum(L==3),sum(L==4)],'descend');
%         %         LL=4*ones(size(idx,1),1);
%         %         LL(L==bi(1))=1;
%         %         LL(L==bi(2))=2;
%         %         LL(L==bi(3))=3;
%         
%         %% roughly nuclei type determination
%         %LL=4*ones(length(idx),1); %
%         
%         %         [mv,idd]=max(C(:,1));  % tumor nuclei: largest size
%         %         LL(idx==idd)=1;        % tumor nuclei
%         %
%         %         mmint=C(:,2);
%         %         mmint(idd)=Inf;
%         %         [mv2,idd2]=min(mmint); % lymphocyte: lower intensity
%         %         LL(idx==idd2)=2;       % lymphocyte
%         
%         %         xx=1:1:numc;
%         %         [mv,idd]=max(C(:,1));  % tumor nuclei: largest size
%         %         LL(idx(:,idd)>0.5)=1;  %% assume 1: tumor cell nuclei
%         %
%         %         mmint=C(:,3);
%         %         mmint(idd)=Inf;
%         %         [mv2,idd2]=min(mmint); % lymphocyte: lower intensity
%         %         LL(idx(:,idd2)>0.8)=2;  %% assume 1: tumor cell nuclei
%         %
%         %         idd3=setdiff(xx,[idd,idd2]);
%         %         LL(idx(:,idd3)>0.8)=3;  %% assume 1: tumor cell nuclei
%         
%         
%         feat_cell=zeros(numel(imgseg),size(ffn_img,2)*fb*2+size(ffn_img,2)*numc+numc);
%         %feat_cell=zeros(numel(imgseg),size(ffn_img,2)*4);
%         inds=1;
%         for ss=1:numel(imgseg)
%             bwn=imread(strcat(seg_path{ip},imgseg(ss).name));
%             temp=strsplit(imgseg(ss).name,'.');
%             img=imread(strcat(image_path{ip},temp{1},'.png'));
%             
%             fnxy2=ffc{ss,1};
%             inde=inds+size(fnxy2,1)-1;
%             
%             
%             %%1) feature bins
%             F00=xu_feature_bins_extraction(ffn_img(inds:inde,:),fb);
%             %F00=xu_feature_bins_extractionII(ffn_img(inds:inde,:));
%             
%             img2=imoverlay(img,bwperim(bwn));
%             figure,imshow(img2);
%             
%             indc=L(inds:inde);
%             bin=1:1:max(L)+1;
%             [bv,bi]=sort(histcounts(indc,bin),'descend'); %%[bv,bi]=sort([sum(indc==1),sum(indc==2),sum(indc==3),sum(indc==4)],'descend');
%             LL=zeros(size(indc,1),1);
%             for k=1:length(bi)
%                 LL(indc==bi(k))=k;
%             end
%             
%             
%             %%2) feature clustering
%             %LL=L(inds:inde);
%             F01=xu_feature_clustering_extractionII(ffn_img(inds:inde,:),LL);
%             %             LL=(numc+1)*ones(size(indc,1),1);    %% for gaussian
%             %             LL(indc==bi(1))=1;
%             %             LL(indc==bi(2))=2;
%             %             LL(indc==bi(3))=3;
%             
%             cnum_str=strsplit(temp{1},'_');
%             if length(cnum_str)>1
%                 rn=str2num(cnum_str{2});
%             else
%                 rn=str2num(cnum_str{1}(24:end));
%             end
%             temp_feat=zeros(1,size(feat_cell,2));
%             temp_feat(1:length([F00,F01]))=[F00,F01];
%             %temp_feat(1:length(F00))=F00;
%             feat_cell(rn,:)=temp_feat;
%             
%             for k=1:max(L(:))
%                 fnxy_temp=fnxy2(LL==k,:);
%                 
%                 switch k
%                     case 1
%                         hold on,plot(fnxy_temp(:,1),fnxy_temp(:,2),'g+');
%                     case 2
%                         hold on,plot(fnxy_temp(:,1),fnxy_temp(:,2),'w+');
%                     case 3
%                         hold on,plot(fnxy_temp(:,1),fnxy_temp(:,2),'b+');
%                     case 4
%                         hold on,plot(fnxy_temp(:,1),fnxy_temp(:,2),'c+');
%                     case 5
%                         hold on,plot(fnxy_temp(:,1),fnxy_temp(:,2),'r+');
%                     case 6
%                         hold on,plot(fnxy_temp(:,1),fnxy_temp(:,2),'y+');
%                     otherwise
%                         disp('impossible!!!!!!!!!!!');
%                 end
%                 
%             end
%             
%             inds=inde+1;
%             
%              saveas(gcf,strcat(debug_path{ip},imgseg(ss).name));
%              close all;
%         end
%         save(strcat(feat_output{ip},cnum_str{1}(1:23),'.mat'),'feat_cell');
        
    end
end