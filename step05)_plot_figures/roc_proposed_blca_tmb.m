%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file script is used for plotting ROC curves for baseline comparison

%  Inputs:
%  .mat files: stores classifier prediction scores, and ground truth labels

%  Output:
%  .fig file: ploting the ROC curves

% Note that: perfcurve() - matlab statistical and machine learning toolbox function
%            line_fewer_markers() - a public available matlab function,
%            download it by google
%            tighgfig() - a public available matlab function, download it
%            by google



% (c) Edited by Hongming Xu,
% Deptment of Quantitative Health Sciences,
% Cleveland Clinic, USA.  December 2017
% If you have any problem feel free to contact me.
% Please address questions or comments to: mxu@ualberta.ca

% Terms of use: You are free to copy,
% distribute, display, and use this work, under the following
% conditions. (1) You must give the original authors credit. (2) You may
% not use or redistribute this work for commercial purposes. (3) You may
% not alter, transform, or build upon this work. (4) For any reuse or
% distribution, you must make clear to others the license terms of this
% work. (5) Any of these conditions can be waived if you get permission
% from the authors.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath(genpath('Y:\projects\xhm_code_repos\matlab_repository\my_codes\'));
matpath='./tcga_blca_proposed/';  %% bladder cancer tmb

load(strcat(matpath,'P_E_TD.mat'));
response_td=labels;
p_e_td=SSC;

load(strcat(matpath,'P_E_APC.mat'))
response_apc=labels;
p_e_apc=SSC;

load(strcat(matpath,'P_E_CN.mat'))
response_cn=labels;
p_e_cn=SSC;

load(strcat(matpath,'p_inceptionv3.mat')); % bladder cancer tmb
response=labels;
p_inceptionv3=SSC;

load(strcat(matpath,'p_resnet50.mat')); % bladder cancer tmb lbp
response_resnet50=labels;
p_resnet50=SSC;

load(strcat('./','p_xception.mat')); % bladder cancer tmb lbp
response_xception=labels;
p_xception=SSC;

% matpath='E:\Hongming\projects\tcga-lung-LUAD\Hongming_codes\4)plot_figures\';  %% lung cancer tmb
% load(strcat(matpath,'score_label_lbp.mat')); % lung cancer tmb
% response_lbp=labels;
% p_gaussian_ss_lbp=SSC;
% 
% load(strcat(matpath,'score_label.mat')); % lung cancer tmb
% response=labels;
% p_gaussian_ss=SSC;







% ROC for g6
[x1_g,y1_g,t1_g,auc1_g]=perfcurve(response_td,p_e_td(:,2),'3');
[x3_g,y3_g,t3_g,auc3_g]=perfcurve(response_cn,p_e_cn(:,2),'3');
[x2_g,y2_g,t2_g,auc2_g]=perfcurve(response_apc,p_e_apc(:,2),'3');

[x4_g,y4_g,t4_g,auc4_g]=perfcurve(response,p_inceptionv3(:,2),'3');
[x5_g,y5_g,t5_g,auc5_g]=perfcurve(response_resnet50,p_resnet50(:,2),'3');
[x6_g,y6_g,t6_g,auc6_g]=perfcurve(response_xception,p_xception(:,2),'3');


figure,
hold on,h1=line_fewer_markers(x1_g,y1_g,7,'p','Color',[0 0.0 0.8],'Spacing', 'curve','markersize',6,'LineWidth',1);
hold on,h2=line_fewer_markers(x3_g,y3_g,7,'x','Color',[0.8 0.0 0],'Spacing', 'curve','markersize',6,'LineWidth',1);
hold on,h3=line_fewer_markers(x2_g,y2_g,7,'*','Color',[0.8 0.5 0],'Spacing', 'curve','markersize',6,'LineWidth',1);

hold on,h4=line_fewer_markers(x4_g,y4_g,7,'d','Color',[0 0.8 0],'Spacing', 'curve','markersize',6,'LineWidth',1);
hold on,h5=line_fewer_markers(x5_g,y5_g,7,'','Color',[0 0.4470 0.7410],'Spacing', 'curve','markersize',6,'LineWidth',1);
hold on,h6=line_fewer_markers(x6_g,y6_g,7,'ks','Spacing', 'curve','markersize',6,'LineWidth',1);

x=0:0.05:1;
y=0:0.05:1;
hold on,h7=plot(x,y,'k--');
legend([h1,h2,h3,h4,h5,h6],{'P-E-TD (AUC=0.683)','P-E-CN (AUC=0.687)','P-E-APC (AUC=0.753)','P-InceptionV3 (AUC=0.713)','P-Resnet50 (AUC=0.747)','P-Xception (AUC=0.752)'});
xlabel('1-specificity'); 
ylabel('sensitivity');
ylim([0 1.0])
title('Low TMB vs High TMB');
grid on;
tightfig();
