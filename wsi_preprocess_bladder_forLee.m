

function [bwTissue,bwHimg2]=wsi_preprocess_bladder_forLee(RGB,thrWhite,thrOverStain,thrNoise)

distob=200;                                                     %% to avoid select boundary patches to againest markers on image borders

%% step 1: select interested regions
gimg=rgb2gray(RGB);
bw=(gimg<=thrWhite);
bwTissue=false(size(bw));
bwTissue(distob:end-distob,distob:end-distob)=bw(distob:end-distob,distob:end-distob);
bwTissue=bwareaopen(bwTissue,thrNoise);                         % obtain tissue pixels


%% step 2--fill small holes in the image--%
bwNoholes=imfill(bwTissue,'holes');
holes=bwNoholes&~bwTissue;
bigholes=bwareaopen(holes,round(thrNoise/5));    % holes greater than half of image patch not filled
smallholes=holes&~bigholes;
bwTissue=bwTissue|smallholes;
%-- end filling small holes -- %

%% step 3: color unmixing to select nuclei channel
bwimgO=(RGB(:,:,3)>thrOverStain);                        % over staining regions

RGB=imresize(RGB,0.5);                                   % for memory issue
[~,Himg,~]=normalizeStaining(RGB);                       
gHimg=rgb2gray(Himg);
gHimg=imresize(gHimg,size(bwTissue));                    % for memory issue
%                            
thresh=multithresh(gHimg(bwimgO),2);
bwHimg=(gHimg<=thresh(1));                              % step 3: obtain nuclei regions

bwHimg=bwHimg&bwimgO;

bwHimg2=false(size(bwHimg));
bwHimg2(distob:end-distob,distob:end-distob)=bwHimg(distob:end-distob,distob:end-distob);
