% matlab script to plot roc curves
% tcga_blca tmb prediction comparison with baseline methods


addpath(genpath('E:/matlab_repository/misc/'));
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


[x0_g,y0_g,t0_g,auc0_g]=perfcurve(response_lbp,p_lbp(:,2),'3','NBoot',1000,'Alpha',0.05); % 95% confidence interval
[x1_g,y1_g,t1_g,auc1_g]=perfcurve(response,p_vgg,'3','NBoot',1000,'Alpha',0.05);
[x2_g,y2_g,t2_g,auc2_g]=perfcurve(response_vgg16_tl,p_vggtl,'3','NBoot',1000,'Alpha',0.05);
[x3_g,y3_g,t3_g,auc3_g]=perfcurve(response,p_xception(:,2),'3','NBoot',1000,'Alpha',0.05);


figure,
hold on,h1=line_fewer_markers(x0_g(:,1),y0_g(:,1),7,'p','Color',[0 0.0 0.8],'Spacing', 'curve','markersize',6,'LineWidth',1);
hold on,h2=line_fewer_markers(x1_g(:,1),y1_g(:,1),7,'d','Color',[0.8 0.0 0],'Spacing', 'curve','markersize',6,'LineWidth',1);
hold on,h3=line_fewer_markers(x2_g(:,1),y2_g(:,1),7,'o','Color',[0.0 0.8 0],'Spacing', 'curve','markersize',6,'LineWidth',1);
hold on,h4=line_fewer_markers(x3_g(:,1),y3_g(:,1),7,'ks','Spacing', 'curve','markersize',6,'LineWidth',1);

x=0:0.05:1;
y=0:0.05:1;
hold on,h7=plot(x,y,'k--');
legend([h1,h2,h3,h4],{'LBP+SVM (AUC=0.623)','SVGG16 (AUC=0.651)','VGG16-TL (AUC=0.707)','Proposed (AUC=0.752)'});
xlabel('1-specificity'); 
ylabel('sensitivity');
ylim([0 1.0])
title('Low TMB vs High TMB');
grid on;
tightfig();

