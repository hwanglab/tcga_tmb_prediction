clear vars;
featPath={'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\2)mutation_prediction\low\','E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\2)mutation_prediction\mid\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\2)mutation_prediction\high\'};

Ftrain=[];
Ltrain=[];

cell_featPath={'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\3)mutation_prediction_cellular_features\low\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\3)mutation_prediction_cellular_features\mid\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\3)mutation_prediction_cellular_features\high\'};

featp={'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\toSunho\low\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\toSunho\mid\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\toSunho\high\'};

% step 1) feature organization
for ip=1:length(featPath)
    fp=featPath{ip};
    cp=cell_featPath{ip};
    fo=featp{ip};
    mats=dir(fullfile(fp,'*.mat'));
    
    for im=1:numel(mats)
        matfilename=strcat(fp,mats(im).name);
        load(matfilename);
        
        coeff=feat_out(:,end)./sum(feat_out(:,end));
        %1) textural features
        feat00=feat_out(:,1:end-1);
        coeff00=repmat(coeff,1,size(feat00,2));
        feat_mean00=sum(feat00.*coeff00,1);
        
%         %2) cellular features
         cmatfilename=strcat(cp,mats(im).name(1:23),'.mat');
         load(cmatfilename);
         coeff01=repmat(coeff,1,size(feat_cell,2));
         feat_mean01=sum(feat_cell.*coeff01,1);
        
        feat_mean=[feat_mean00,feat_mean01(1:70)];
        Ftrain=[Ftrain;feat_mean];
        
        if ip==1
            Llow=ones(size(feat_mean,1),1);
            Ltrain=[Ltrain;Llow];
        elseif ip==2
            Lmid=2*ones(size(feat_mean,1),1);
            Ltrain=[Ltrain;Lmid];
        else
            Lhigh=3*ones(size(feat_mean,1),1);
            Ltrain=[Ltrain;Lhigh];
        end
        
        save(strcat(fo,mats(im).name),'feat_mean');
    end
    
end