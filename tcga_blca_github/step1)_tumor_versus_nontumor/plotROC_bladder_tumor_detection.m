% plot ROC curve for bladder cancer detection

load('bladder_score_label.mat'); % bladder cancer
response=testingResponse;
p_gaussian_ss=foldScores;

[x6_g,y6_g,t6_g,auc6_g]=perfcurve(response,p_gaussian_ss(:,2),'1');


figure,
hold on,h6=line_fewer_markers(x6_g,y6_g,7,'bs','Spacing', 'curve','markersize',8,'LineWidth',1.5);

x=0:0.05:1;
y=0:0.05:1;
hold on,h7=plot(x,y,'k--');
legend([h6,h7],{'LBP+SVM (AUC=0.98)','Reference Line'});
xlabel('1-Specificity'); 
ylabel('Sensitivity');
ylim([0 1.0])
title('Tumor versus Non-Tumor');
grid on;


tightfig();