addpath(genpath('E:\matlab_repository\misc\'));
matpath='E:\Hongming\projects\tcga-bladder-mutationburden\Hongming_codes\step7)_plot_figures\';  %% bladder cancer tmb

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

load(strcat(matpath,'p_xception.mat')); % bladder cancer tmb lbp
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
hold on,h6=line_fewer_markers(x6_g,y6_g,7,'ks','Spacing', 'curve','markersize',6,'LineWidth',1.5);

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
% % ROC for g7
% [x7_fa_c,y7_fa_c,t7_fa_c,auc7_fa_c]=perfcurve(label,fa_cubic_ss(:,2),'7');
% [x7_fa_g,y7_fa_g,t7_fa_g,auc7_fa_g]=perfcurve(label,fa_gaussian_ss(:,2),'7');
% 
% 
% [x7_hhg_g,y7_hhg_g,t7_hhg_g,auc7_hhg_g]=perfcurve(label,hhg_gaussian_ss(:,2),'7');
% 
% [x7_c,y7_c,t7_c,auc7_c]=perfcurve(label,p_cubic_ss(:,2),'7');
% [x7_g,y7_g,t7_g,auc7_g]=perfcurve(label,p_gaussian_ss(:,2),'7');
% 
% 
% figure,
% hold on,h1=line_fewer_markers(x7_fa_c,y7_fa_c,7,'m','Spacing', 'curve','markersize',8,'LineWidth',1.5);
% hold on,h2=line_fewer_markers(x7_fa_g,y7_fa_g,7,'cp','Spacing', 'curve','markersize',8,'LineWidth',1.5);
% 
% hold on,h3=line_fewer_markers(x7_hhg_c,y7_hhg_c,7,'gd','Spacing', 'curve','markersize',8,'LineWidth',1.5);
% hold on,h4=line_fewer_markers(x7_hhg_g,y7_hhg_g,7,'k*','Spacing', 'curve','markersize',8,'LineWidth',1.5);
% 
% hold on,h5=line_fewer_markers(x7_c,y7_c,7,'rx','Spacing', 'curve','markersize',8,'LineWidth',1.5);
% hold on,h6=line_fewer_markers(x7_g,y7_g,7,'bs','Spacing', 'curve','markersize',8,'LineWidth',1.5);
% 
% x=0:0.05:1;
% y=0:0.05:1;
% hold on,h7=plot(x,y,'k--');
% legend([h1,h2,h3,h4,h5,h6],{'FA-Polynomial','FA-Gaussian','HHG-Polynomial','HHG-Gaussian','CSLBP-Polynomial','CSLBP-Gaussian'});
% xlabel('False positive rate'); 
% ylabel('True positive rate');
% ylim([0 1.03])
% title('intermediate risk VS low and high risks');
% grid on;
% 
% 
% 
% % ROC for g8
% [x8_fa_c,y8_fa_c,t8_fa_c,auc8_fa_c]=perfcurve(label,fa_cubic_ss(:,3),'8');
% [x8_fa_g,y8_fa_g,t8_fa_g,auc8_fa_g]=perfcurve(label,fa_gaussian_ss(:,3),'8');
% 
% [x8_hhg_c,y8_hhg_c,t8_hhg_c,auc8_hhg_c]=perfcurve(label,hhg_cubic_ss(:,3),'8');
% [x8_hhg_g,y8_hhg_g,t8_hhg_g,auc8_hhg_g]=perfcurve(label,hhg_gaussian_ss(:,3),'8');
% 
% [x8_c,y8_c,t8_c,auc8_c]=perfcurve(label,p_cubic_ss(:,3),'8');
% [x8_g,y8_g,t8_g,auc8_g]=perfcurve(label,p_gaussian_ss(:,3),'8');
% 
% 
% figure,
% hold on,h1=line_fewer_markers(x8_fa_c,y8_fa_c,7,'m','Spacing', 'curve','markersize',8,'LineWidth',1.5);
% hold on,h2=line_fewer_markers(x8_fa_g,y8_fa_g,7,'cp','Spacing', 'curve','markersize',8,'LineWidth',1.5);
% 
% hold on,h3=line_fewer_markers(x8_hhg_c,y8_hhg_c,7,'gd','Spacing', 'curve','markersize',8,'LineWidth',1.5);
% hold on,h4=line_fewer_markers(x8_hhg_g,y8_hhg_g,7,'k*','Spacing', 'curve','markersize',8,'LineWidth',1.5);
% 
% hold on,h5=line_fewer_markers(x8_c,y8_c,7,'rx','Spacing', 'curve','markersize',8,'LineWidth',1.5);
% hold on,h6=line_fewer_markers(x8_g,y8_g,7,'bs','Spacing', 'curve','markersize',8,'LineWidth',1.5);
% 
% x=0:0.05:1;
% y=0:0.05:1;
% hold on,h7=plot(x,y,'k--');
% legend([h1,h2,h3,h4,h5,h6],{'FA-Polynomial','FA-Gaussian','HHG-Polynomial','HHG-Gaussian','CSLBP-Polynomial','CSLBP-Gaussian'});
% xlabel('False positive rate'); 
% ylabel('True positive rate');
% ylim([0 1.03])
% title('high risk VS low and intermediate risks');
% grid on;

% 
% auc_fa_c=(auc6_fa_c+auc7_fa_c+auc8_fa_c)/3;
% auc_fa_g=(auc6_fa_g+auc7_fa_g+auc8_fa_g)/3;
% 
% auc_hhg_c=(auc6_hhg_c+auc7_hhg_c+auc8_hhg_c)/3;
% auc_hhg_g=(auc6_hhg_g+auc7_hhg_g+auc8_hhg_g)/3;
% 
% auc_cslbp_c=(auc6_c+auc7_c+auc8_c)/3;
% auc_cslbp_g=(auc6_g+auc7_g+auc8_g)/3;