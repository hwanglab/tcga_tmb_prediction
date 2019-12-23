%% main function for selecting patches for deep learning training

close all;clc;

imagePath={'E:\blca_mutationBurden\blca_wsi\'};

indices_output='Y:\DL_TMB_prediction\trail01\';

%% (1) build a cross validation index
load('E:\Hongming\projects\tcga-bladder-mutationburden\Hongming_codes\blca_MutBurdens.mat');  % TCGA TMB values
wsi_label=[];
for ip=1:length(imagePath)
    
    imgPath=imagePath{ip};
    
    imgs=dir(fullfile(imgPath,'*.svs'));
    
    pId=strings(numel(imgs),1);
    
    for im=1:numel(imgs)
        patientID=imgs(im).name(1:23);
        temp=strcmp(blca_MutBurdens(:,1),patientID(1:12));
        
        pId(im,1)=patientID(1:12);
        
        catg=blca_MutBurdens(temp==1,3);
        if strcmp(catg,'Low')
            wsi_label=[wsi_label;10];
        elseif strcmp(catg,'Mid')
            wsi_label=[wsi_label;20];
        elseif strcmp(catg,'High')
            wsi_label=[wsi_label;30];
        else
            disp('impossible!!!!');
        end
    end
end
rng('default');
rng(1);
KFolds=3;
Indices=crossvalind('Kfold',wsi_label,KFolds);

cvIndices=zeros(numel(blca_MutBurdens(:,1)),1);
for jp=1:numel(blca_MutBurdens(:,1))
    temp=strcmp(blca_MutBurdens{jp,1},pId);
    if sum(temp)==1
        cvIndices(jp,1)=Indices(temp==1);
    end
end


% for the other 24 images
imagePath2={'E:\blca_mutationBurden\tumor_detection\'};
wsi_label=[];
for ip=1:length(imagePath2)
    
    imgPath=imagePath2{ip};
    
    imgs=dir(fullfile(imgPath,'*.svs'));
    
    pId=strings(numel(imgs),1);
    
    for im=1:numel(imgs)
        patientID=imgs(im).name(1:23);
        temp=strcmp(blca_MutBurdens(:,1),patientID(1:12));
        
        pId(im,1)=patientID(1:12);
        
        catg=blca_MutBurdens(temp==1,3);
        if strcmp(catg,'Low')
            wsi_label=[wsi_label;10];
        elseif strcmp(catg,'Mid')
            wsi_label=[wsi_label;20];
        elseif strcmp(catg,'High')
            wsi_label=[wsi_label;30];
        else
            disp('impossible!!!!');
        end
    end
end
rng('default');
rng(1);
KFolds=3;
Indices=crossvalind('Kfold',wsi_label,KFolds);

for jp=1:numel(blca_MutBurdens(:,1))
    temp=strcmp(blca_MutBurdens{jp,1},pId);
    if sum(temp)==1
        cvIndices(jp,1)=Indices(temp==1);
    end
end

save(strcat(indices_output,'cvInd.mat'),'cvIndices')


