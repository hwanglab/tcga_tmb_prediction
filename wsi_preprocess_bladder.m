

function [bwTissue,bwHimg2]=wsi_preprocess_bladder(RGB,thrWhite,thrOverStain)

distob=220;                                                     %% to avoid select boundary patches to againest markers on image borders

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
bigholes=bwareaopen(holes,round(thrNoise/10));    % holes greater than half of image patch not filled
smallholes=holes&~bigholes;
bwTissue=bwTissue|smallholes;
%-- end filling small holes -- %

%% step 3: color unmixing to select nuclei channel
bwimgO=(RGB(:,:,3)>thrOverStain);                        % over staining regions

% RGB=imresize(RGB,0.5);                                   % for memory issue
% [~,Himg,~]=normalizeStaining(RGB);                       
% gHimg=rgb2gray(Himg);
% gHimg=imresize(gHimg,size(bwTissue));                    % for memory issue
%                            
% thresh=multithresh(gHimg(bwimgO),2);
% bwHimg=(gHimg<=thresh(1));                              % step 3: obtain nuclei regions
t_nuclei=0.7*255;
bwHimg=(RGB(:,:,1)<t_nuclei);
bwHimg=bwHimg&bwimgO;

bwHimg2=false(size(bwHimg));
bwHimg2(distob:end-distob,distob:end-distob)=bwHimg(distob:end-distob,distob:end-distob);
