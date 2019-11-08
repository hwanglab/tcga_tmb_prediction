clear vars;
featPath={'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\2)mutation_prediction\low\','E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\2)mutation_prediction\mid\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\2)mutation_prediction\high\'};

featmatoutput='E:\Hongming\projects\tcga-bladder-mutationburden\Hongming_codes\step3)_wis_feature_organization\';
Ftrain=[];
Ltrain=[];


% v-1: feature organization for clustering-based methods
for ip=1:length(featPath)
    fp=featPath{ip};
    mats=dir(fullfile(fp,'*.mat'));
    
    for im=1:numel(mats)
        matfilename=strcat(fp,mats(im).name);
        load(matfilename);
        feat00=feat_out(:,1:end-1);
        coeff=feat_out(:,end)./sum(feat_out(:,end));
        coeff00=repmat(coeff,1,size(feat00,2));
        
        feat_mean=sum(feat00.*coeff00,1);
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

% % v-2: feature organization for averaging methods
% for ip=1:length(featPath)
%     fp=featPath{ip};
%     mats=dir(fullfile(fp,'*.mat'));
%     
%     for im=1:numel(mats)
%         matfilename=strcat(fp,mats(im).name);
%         load(matfilename)
%         
%         feat_mean=mean(feat_np);
%         %feat_std=std(feat_np);
%         %feat_kur=kurtosis(feat_tumor);
%         %feat_skew=skewness(feat_tumor);
%         Ftrain=[Ftrain;feat_mean];
%         
%         if ip==1
%             Llow=ones(size(feat_mean,1),1);
%             Ltrain=[Ltrain;Llow];
%         elseif ip==2
%             Lmid=2*ones(size(feat_mean,1),1);
%             Ltrain=[Ltrain;Lmid];
%         else
%             Lhigh=3*ones(size(feat_mean,1),1);
%             Ltrain=[Ltrain;Lhigh];
%         end
%     end
% end


FT=table(Ftrain,Ltrain);
FT.Properties.VariableNames={'features','classes'};

saveName=strcat(featmatoutput,'FT.mat');
save(saveName,'FT')