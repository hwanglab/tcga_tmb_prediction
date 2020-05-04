

function [bwTissue]=wsi_preprocess_foreground(RGB,thrWhite)

distob=200;                                                     %% to avoid select boundary patches to againest markers on image borders

%% step 1: select interested regions
gimg=rgb2gray(RGB);
bw=(gimg<=thrWhite);
bwTissue=false(size(bw));
bwTissue(distob:end-distob,distob:end-distob)=bw(distob:end-distob,distob:end-distob);

CC=bwconncomp(bwTissue);
numPixels = cellfun(@numel,CC.PixelIdxList);
thrNoise=round(max(numPixels)/5);                               % less than 1/5 of the maximum region is considered as noise
bwTissue=bwareaopen(bwTissue,thrNoise);                         % obtain tissue pixels


%% step 2--fill small holes in the image--%
bwNoholes=imfill(bwTissue,'holes');
holes=bwNoholes&~bwTissue;
bigholes=bwareaopen(holes,round(thrNoise/5));    % holes greater than half of image patch not filled
smallholes=holes&~bigholes;
bwTissue=bwTissue|smallholes;
%-- end filling small holes -- %

