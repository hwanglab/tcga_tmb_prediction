% for generating the tissue mask for whole slide .svs image
% input: RGB->.svs slide image
%        thrWhite-> threshold, i.e., pixels above this value are considered
%        as white background pixels
%        tilesize-> threshold, i.e., holes above this value are considered
%        as big holes, not filled by morphological operations
% output: bwTissue-> binary mask with foreground as tissue regions
% author: Hongming Xu, Cleveland Clinic
% feel free to use it, but need give original author credits

function [bwTissue]=wsi_preprocess_tissue(RGB,thrWhite,tilesize)

distob=200;                                                     %% to avoid select boundary patches to againest markers on image borders

%% step 1: select interested regions
gimg=rgb2gray(RGB);
bw=(gimg<=thrWhite);
bwTissue=false(size(bw));
bwTissue(distob:end-distob,distob:end-distob)=bw(distob:end-distob,distob:end-distob);

CC=bwconncomp(bwTissue);
numPixels = cellfun(@numel,CC.PixelIdxList);
thrNoise=round(max(numPixels)/2);                               % less than half of the maximum region is considered as noise
bwTissue=imclose(bwTissue,strel('disk',2));
bwTissue=bwareaopen(bwTissue,thrNoise);                         % obtain tissue pixels


%% step 2--fill small holes in the image--%
bwNoholes=imfill(bwTissue,'holes');
holes=bwNoholes&~bwTissue;
bigholes=bwareaopen(holes,tilesize);    % holes greater than half of image patch not filled
smallholes=holes&~bigholes;
bwTissue=bwTissue|smallholes;
%-- end filling small holes -- %

