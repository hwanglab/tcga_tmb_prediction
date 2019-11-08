% Main file to run different color normalization techniques

clear
clc
close all
addpath(genpath('E:/matlab_repository/toolboxes/spams-matlab-v2.6-2017-02-27/spams-matlab-v2.6/'));
addpath(genpath('E:/Hongming/resources/CodeRelease_ColorNormalization-master/SNMF stain separation and color normalization/'));

%% Define source and target images (Add target and source images to the folder "images")
sourcepath='E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\tissue_images_normalized\';
targetpath='E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\tissue_images_normalized\TCGA-E2-A14V-01Z-00-DX1.png';
%source=imread('images/source1.png');
target=imread(targetpath);
temp=strsplit(targetpath,'\');

nstains=2;
lambda=0.1;  % Use smaller values of the lambda (0.01-0.1) for better reconstruction. however, if the normalized image seems not fine, increase or decrease the value accordingly.

%% Our Method (The techniques is published in ISBI 2015 under the title "STRUCTURE-PRESERVED COLOR NORMALIZATION FOR HISTOLOGICAL IMAGES")
% For queries, contact: abhishek.vahadane@gmail.com, vahadane@iitg.ernet.in
% Source and target stain separation and storage of factors
tic
[Wi, Hi,Hiv]=stainsep(target,nstains,lambda);
% save('target.mat','Wi','Hi','Hiv')

list=dir(strcat(sourcepath,'*.png'));

for k=1:length(list)
    filename=list(k).name;
    if ~strcmp(filename,temp{end})
        source=imread(strcat(sourcepath,filename));
        [Wis, His,Hivs]=stainsep(source,nstains,lambda);
        % save('source.mat','Wis','His','Hivs')
        
        
        % Color normlization
        % addpath(genpath('Our Method'))
        [our]=SCN(source,Hi,Wi,His);
        time=toc;
        
        % Write image to a folder
        %imwrite(our,'images/norm.png')
        %% Visuals
%         figure;
%         subplot(131);imshow(source);xlabel('source')
%         subplot(132);imshow(target);xlabel('target')
%         subplot(133);imshow(our);xlabel('normalized source')
        imwrite(our,strcat(sourcepath,filename(1:23),'.png'));
    else
        disp('target image!!!!');
        our=imread(strcat(sourcepath,filename));
        imwrite(our,strcat(sourcepath,filename(1:23),'.png'));
    end
    
end