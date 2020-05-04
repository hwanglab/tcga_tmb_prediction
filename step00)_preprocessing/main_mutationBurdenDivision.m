
%1) bladder cancer patients tcga TMB table
% % read and process excel spreadsheet file
% sheet=2;
% [~,~,raw_patientID]=xlsread('Table_S1.2017_08_05.xlsx',sheet,'A1:A413');
% [~,~,raw_MB]=xlsread('Table_S1.2017_08_05.xlsx',sheet,'DK1:DK413');
% MBm=cell2mat(raw_MB(2:end));
% 
% LI=33;
% IH=66;
% tLI=prctile(MBm,LI);
% tIH=prctile(MBm,IH);
% 
% blca_MutBurdens=cell(length(raw_patientID)-1,3);
% blca_MutBurdens(:,1)=raw_patientID(2:end);
% blca_MutBurdens(:,2)=raw_MB(2:end);
% blca_MutBurdens(MBm<=tLI,3)={'Low'};
% blca_MutBurdens(MBm>tLI & MBm<=tIH,3)={'Mid'};
% blca_MutBurdens(MBm>tIH,3)={'High'};
% 
% save('blca_MutBurdens.mat','blca_MutBurdens');

%2) bladder cancer patients Yunku provided TMB values
sheet=1;
[~,~,raw_patientID]=xlsread('blca_TMB.xlsx',sheet,'A1:A413');
[~,~,raw_MB]=xlsread('blca_TMB.xlsx',sheet,'B1:B413');
MBm=cell2mat(raw_MB(2:end));

LI=33;
IH=66;
tLI=prctile(MBm,LI);
tIH=prctile(MBm,IH);

blca_MutBurdens=cell(length(raw_patientID),3);
temp=cell2mat(raw_patientID(2:end));
temp2=temp(:,1:12);
temp3=mat2cell(temp2,ones(size(temp2,1),1),12);
blca_MutBurdens(:,1)=temp3;
blca_MutBurdens(:,2)=raw_MB(2:end);
blca_MutBurdens(MBm<=tLI,3)={'Low'};
blca_MutBurdens(MBm>tLI & MBm<=tIH,3)={'Mid'};
blca_MutBurdens(MBm>tIH,3)={'High'};

save('blca_MutBurdens_Yunku.mat','blca_MutBurdens');

% low: 137 patients, intermediate: 135 patients, high: 140 patients