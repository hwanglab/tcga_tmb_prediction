% evaluation for tcga_luad tmb prediction
% result reported in the paper
% author: Hongming Xu, CCF, 2019
% questions: mxu@ualberta.ca


function VGG16_BLCA

disp('loading patient ID and ground truth TMB levels!')

load('./patID253.mat')

dl=1; % dl-method: svgg16
tl=1; % tl-method: vgg16 transfer learning

if dl==1
    disp('loading automatic predictions from SVGG model!')
    load('./svgg16/vggresult.mat');
    
    label=[];
    for i=1:length(patID(:,1))
        if strcmp(patID{i,2},'Low')
            label=[label,1];
        elseif strcmp(patID{i,2},'High')
            label=[label,3];
        else
            disp('impossible!!!!');
        end
    end
    
    pp=3*ones(1,length(label));
    pp(preds>0.46)=1;
    
    ACC=sum(label==pp)/numel(label);
    SPE=sum(label(label==1)==pp(label==1))/(sum(label==1));
    SNE=sum(label(label==3)==pp(label==3))/(sum(label==3));
    
    
    % only used for KM-curve plotting
    pred_bladder=cell(length(pp),3);
    pred_bladder(:,1)=patID(:,1);
    pred_bladder(pp==1,2)={'Low'};
    pred_bladder(pp==3,2)={'High'};
    pred_bladder(label==1,3)={'Low'};
    pred_bladder(label==3,3)={'High'};
end

if tl==1
    load('./vgg16_tl/tl_vgg16.mat')
    
    labels=[];
    for i=1:length(pID)
        temp_p=pID(i,:);
        ind=find(strcmp(temp_p,patID(:,1)));
        if strcmp(patID{ind,2},'Low')
            labels=[labels,1];
        elseif strcmp(patID{ind,2},'High')
            labels=[labels,3];
        else
            disp('impossible!!!!');
        end
    end
    
    pp=3*ones(1,length(labels));
    pp(preds>0.5)=1;
    
    ACC=sum(labels==pp)/numel(labels);
    SPE=sum(labels(labels==1)==pp(labels==1))/(sum(labels==1));
    SNE=sum(labels(labels==3)==pp(labels==3))/(sum(labels==3));
    
    %save('vgg16_tl_roc.mat','preds','labels');
end

