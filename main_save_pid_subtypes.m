%1) bladder cancer patients tcga TMB table
% read and process excel spreadsheet file
sheet=2;
[~,~,raw_patientID]=xlsread('Table_S1.2017_08_05.xlsx',sheet,'A1:A413');
[~,~,raw_MB]=xlsread('Table_S1.2017_08_05.xlsx',sheet,'AV1:AV413');

blca_Subtypes=cell(length(raw_patientID)-1,2);
blca_Subtypes(:,1)=raw_patientID(2:end);
blca_Subtypes(:,2)=raw_MB(2:end);
% 
save('blca_Subtypes.mat','blca_Subtypes');