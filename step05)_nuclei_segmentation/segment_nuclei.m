function segment_nuclei

% mask_path={'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_masks\low\',...
%     'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_masks\mid\',...
%     'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_masks\high\'};
% 
% image_path={'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_images\low\',...
%     'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_images\mid\',...
%     'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_images\high\'};
% 
% debug_path={'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_debugs\low\',...
%     'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_debugs\int\',...
%     'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_debugs\high\'};
% 
% seg_path={'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_segs\low\',...
%     'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_segs\int\',...
%     'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_segs\high\'};


mask_path={'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_masksII\low\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_masksII\mid\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_masksII\high\'};

image_path={'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_imagesII\low\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_imagesII\mid\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_imagesII\high\'};

debug_path={'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_debugs\low\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_debugs\int\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_debugs\high\'};

seg_path={'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_segsII\low\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_segsII\mid\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\mutation_analysis_segsII\high\'};

thr=0.8*255;
thrs=0.2*255;
siz_nuc=40;
size_nuc_max=2500;

for mp=1:length(mask_path)
    list=dir(strcat(mask_path{mp},'*.png'));
    for k=1:length(list)
        mask_name=list(k).name;
        mask=imread(strcat(mask_path{mp},mask_name));
        
        img=imread(strcat(image_path{mp},mask_name(1:end-9)));
        bwn=(mask(:,:,1)<thr);
        bwn=bwareaopen(bwn,siz_nuc);
        bwn=imdilate(bwn,strel('disk',1));
        
        bwNoholes=imfill(bwn,'holes');
        holes=bwNoholes&~bwn;
        bigholes=bwareaopen(holes,siz_nuc);    % holes greater than half of image patch not filled
        smallholes=holes&~bigholes;
        bwn=bwn|smallholes;
        
        
        bws=(mask(:,:,3)>thrs);
        temp=mask(:,:,2)>mask(:,:,2);
        bws=bws&(~temp);
        %bws=imopen(bws,strel('disk',1));
        %bws=imfill(bws,'holes');
        bws=bwareaopen(bws,round(siz_nuc/10));
        ssn=bwn&bws;
        %     figure,imshow(bwn);
        %     [r,c]=find(ssn);
        %     hold on,plot(c,r,'r.');
        
        D=bwdist(~bwn);
        D=-D;
        D=imimposemin(D,ssn);
        L=watershed(D);
        bwn(L==0)=0;
        bwn=xor(bwareaopen(bwn,siz_nuc),bwareaopen(bwn,size_nuc_max));
        bwn=imfill(bwn,'holes');
        
        imwrite(bwn,strcat(seg_path{mp},mask_name));
        
        img2=imoverlay(img,bwperim(bwn));
        figure,imshow(img2);
        %hold on,plot(c,r,'r.');
        saveas(gcf,strcat(debug_path{mp},mask_name));
        close all;
    end
end