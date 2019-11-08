clear vars;
addpath(genpath('E:\Hongming\projects\tcga-bladder-mutationburden\Hongming_codes\step5)_nuclei_segmentation\'));
featPath={'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\6)mutation_prediction\low\','E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\6)mutation_prediction\mid\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\6)mutation_prediction\high\'};

featmatoutput='E:\Hongming\projects\tcga-bladder-mutationburden\Hongming_codes\step3)_wis_feature_organization\';
Ftrain=[];
Ltrain=[];
fb=10;

% v-2: feature organization for averaging methods
for ip=1:length(featPath)
    fp=featPath{ip};
    mats=dir(fullfile(fp,'*.mat'));
    
    for im=1:numel(mats)
        matfilename=strcat(fp,mats(im).name);
        load(matfilename)
        
        %feat_mean=xu_feature_bins_extraction(feat_out(1:112),fb);
        feat_mean=mean(feat_out(:,1:112));
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
    end
end


FT=table(Ftrain,Ltrain);
FT.Properties.VariableNames={'features','classes'};

saveName=strcat(featmatoutput,'FT.mat');
save(saveName,'FT')