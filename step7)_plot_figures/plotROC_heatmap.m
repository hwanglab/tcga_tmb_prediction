addpath(genpath('E:\matlab_repository\misc\'));
matpath='E:\Hongming\projects\tcga-bladder-mutationburden\Hongming_codes\step7)_plot_figures\';  %% bladder cancer tmb

load(strcat(matpath,'bladder_heatmap.mat')); % bladder cancer tmb
response=labels;
p_xception=SSC;

[x2_g,y2_g,t2_g,auc2_g]=perfcurve(response,p_xception(:,2),'3');


figure,
hold on,h1=line_fewer_markers(x2_g,y2_g,7,'ks','Spacing', 'curve','markersize',6,'LineWidth',2);

x=0:0.05:1;
y=0:0.05:1;
hold on,h7=plot(x,y,'k--');
legend(h1,{'Proposed (AUC=0.747)'});
xlabel('1-specificity'); 
ylabel('sensitivity');
ylim([0 1.0])
title('Low TMB vs High TMB');
grid on;
tightfig();