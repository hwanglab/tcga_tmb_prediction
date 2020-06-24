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


addpath(addpath(genpath('..\utility_funcs\')))
matpath='./tcga_blca_baselines/';  

%1) proposed tmb prediction
load(strcat('./','p_xception.mat')); 
response=labels;
p_xception=SSC;

%2) lbp baseline comparison
load(strcat(matpath,'lbp.mat')); % bladder cancer tmb lbp
response_lbp=labels;
p_lbp=SSC;

%3) svgg16 baseline comparison
load(strcat(matpath,'vggresult.mat'))
p_vgg=1-preds;

%4) vgg16-tl baseline comparison
load(strcat(matpath,'vgg16_tl_roc.mat'))
response_vgg16_tl=labels;
p_vggtl=1-preds;

%5) mil baseline comparison
load(strcat(matpath,'mil.mat'))
response_mil=labels_mil;
p_mil=preds_mil;

%6) resnet18 baseline comparison
load(strcat(matpath,'resnet18.mat'))
response_resnet18=label;
p_resnet18=preds;

[x0_g,y0_g,t0_g,auc0_g]=perfcurve(response_lbp,p_lbp(:,2),'3','NBoot',1000,'Alpha',0.05); % 95% confidence interval
[x1_g,y1_g,t1_g,auc1_g]=perfcurve(response,p_vgg,'3','NBoot',1000,'Alpha',0.05);
[x2_g,y2_g,t2_g,auc2_g]=perfcurve(response_vgg16_tl,p_vggtl,'3','NBoot',1000,'Alpha',0.05);
[x3_g,y3_g,t3_g,auc3_g]=perfcurve(response_mil,p_mil,'3','NBoot',1000,'Alpha',0.05);
[x4_g,y4_g,t4_g,auc4_g]=perfcurve(response_resnet18,p_resnet18,'3','NBoot',1000,'Alpha',0.05);

[x5_g,y5_g,t5_g,auc5_g]=perfcurve(response,p_xception(:,2),'3','NBoot',1000,'Alpha',0.05);


figure,
hold on,h1=line_fewer_markers(x0_g(:,1),y0_g(:,1),7,'p','Color',[0 0.0 0.8],'Spacing', 'curve','markersize',6,'LineWidth',1);
hold on,h2=line_fewer_markers(x1_g(:,1),y1_g(:,1),7,'d','Color',[0.8 0.0 0],'Spacing', 'curve','markersize',6,'LineWidth',1);
hold on,h3=line_fewer_markers(x2_g(:,1),y2_g(:,1),7,'o','Color',[0.0 0.8 0],'Spacing', 'curve','markersize',6,'LineWidth',1);
hold on,h4=line_fewer_markers(x3_g(:,1),y3_g(:,1),7,'o','Color',[0 0.4470 0.7410],'Spacing', 'curve','markersize',6,'LineWidth',1);
hold on,h5=line_fewer_markers(x4_g(:,1),y4_g(:,1),7,'*','Color',[0.8 0.5 0],'Spacing', 'curve','markersize',6,'LineWidth',1);

hold on,h6=line_fewer_markers(x5_g(:,1),y5_g(:,1),7,'ks','Spacing', 'curve','markersize',6,'LineWidth',1);

x=0:0.05:1;
y=0:0.05:1;
hold on,h7=plot(x,y,'k--');
%legend([h2,h3,h4],{'VGG16-DL (AUC=0.651)','VGG16-TL (AUC=0.707)','Proposed (AUC=0.752)'});
legend([h1,h2,h3,h4,h5,h6],{'LBP+SVM (AUC=0.623)','Designed CNN (AUC=0.651)','VGG16-TL2 (AUC=0.707)','MIL (AUC=0.647)','Resnet18 (AUC=0.701)','Proposed (AUC=0.752)'});
xlabel('1-specificity'); 
ylabel('sensitivity');
ylim([0 1.0])
title('Low TMB vs High TMB');
grid on;
tightfig();

