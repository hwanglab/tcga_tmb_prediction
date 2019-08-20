%% Author: Hongming Xu, Postdoc, Hwang Lab 2018
% function usage: organize features for classification by SVM classifer APP
% the features are save into .mat file
% V0.0: test and use for TCGA bladder pathology slides

clear vars;
featureoutput={'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\1)tumor_prediction\tumor\',...
    'E:\Hongming\projects\tcga-bladder-mutationburden\featureoutput\1)tumor_prediction\nontumor\'};

Ftrain=[];
Ltrain=[];

featmatoutput='E:\Hongming\projects\tcga-bladder-mutationburden\Hongming_codes\step1)_tumor_versus_nontumor\';

temp2=[];
for ip=1:length(featureoutput)
    featpath=featureoutput{ip};
    mats=dir(fullfile(featpath,'*.mat'));
    
    for im=1:numel(mats)
        matfilename=strcat(featpath,mats(im).name);
        %fprintf('load filename=%s\n%d',mats(im).name,im);
        load(matfilename);
        Ftrain=[Ftrain;feat];
        temp2=[temp2,size(feat,1)];
        
        if ip==1 % tumor
            Ltumor=ones(size(feat,1),1);
            Ltrain=[Ltrain;Ltumor];
        else
            Lnontumor=zeros(size(feat,1),1);
            Ltrain=[Ltrain;Lnontumor];
        end
    end
end

FT=table(Ftrain,Ltrain);
FT.Properties.VariableNames={'features','classes'};

%saveName=strcat(featmatoutput,'FT.mat');
%save(saveName,'FT')

