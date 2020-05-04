% evaluation for tcga_blca tmb prediction
% result reported in the paper
% author: Hongming Xu, CCF, 2019
% questions: mxu@ualberta.ca


function Resnet18_MIL_BLCA

disp('loading patient ID and ground truth TMB levels!')

load('./patID253.mat')

resnet18=0; % resent18 method
mil=1;      % multiple instance learning method

if resnet18==1
    disp('loading automatic predictions from resnet18 model!')
    pred_path='./resnet18/';
    
    label=[];
    pids=[];
    preds=[];
    for i=1:length(patID(:,1))
        pp_id=patID{i,1};
        pred_file=dir(fullfile(pred_path,strcat(pp,'*.xlsx')));
        [num,txt,raw] = xlsread(strcat(pred_path,pred_file.name)); 
        preds=[preds,mean(num(:,3))];
        pids=[pids;pp_id]
        if strcmp(patID{i,2},'Low')
            label=[label,1];
        elseif strcmp(patID{i,2},'High')
            label=[label,3];
        else
            disp('impossible!!!!');
        end
    end
    
    pp=ones(1,length(label));
    pp(preds>0.5)=3;
    
    ACC=sum(label==pp)/numel(label);
    SPE=sum(label(label==1)==pp(label==1))/(sum(label==1));
    SNE=sum(label(label==3)==pp(label==3))/(sum(label==3));
    
    save(strcat(pred_path,'resnet18.mat'),'pids','label','preds')
    
    % only used for KM-curve plotting
    %pred_bladder=cell(length(pp),3);
    %pred_bladder(:,1)=patID(:,1);
    %pred_bladder(pp==1,2)={'Low'};
    %pred_bladder(pp==3,2)={'High'};
    %pred_bladder(label==1,3)={'Low'};
    %pred_bladder(label==3,3)={'High'};
end

if mil==1
    pred_path='./mil/';
    [num,txt,raw]=xlsread(strcat(pred_path,'tcga_blca_all.xlsx')); 
    
    labels_mil=[];
    pids_mil=[];
    preds_mil=[];
    for i=1:length(patID(:,1))
        pp_id=patID{i,1};
        ind=find(strcmp(raw(:,2),pp_id));
        preds_mil=[preds_mil,raw{ind,3}];
        pids_mil=[pids_mil;pp_id];
        if strcmp(patID{i,2},'Low')
            labels_mil=[labels_mil,1];
        elseif strcmp(patID{i,2},'High')
            labels_mil=[labels_mil,3];
        else
            disp('impossible!!!!');
        end
    end
    
    pp=ones(1,length(labels_mil));
    pp(preds_mil>0.7)=3;
    
    ACC=sum(labels_mil==pp)/numel(labels_mil);
    SPE=sum(labels_mil(labels_mil==1)==pp(labels_mil==1))/(sum(labels_mil==1));
    SNE=sum(labels_mil(labels_mil==3)==pp(labels_mil==3))/(sum(labels_mil==3));
    
    save(strcat(pred_path,'mil.mat'),'pids_mil','labels_mil','preds_mil')
end

