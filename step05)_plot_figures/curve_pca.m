x=[20,40,60,80,100,120];
y_inceptionv3=[70.75,59.68,64.03,62.45,65.22,60.08];
y_resnet50=[71.94,62.45,62.85,62.06,64.03,63.24];
y_xception=[64.23,64.43,68.77,68.38,73.12,68.77];

figure,
plot(x,y_inceptionv3,'r-s',x,y_resnet50,'g-*',x,y_xception,'b-o','LineWidth',2,'MarkerSize',10)
xticks([20,40,60,80,100,120])
xlabel('number of PCA components');
ylabel('classification accuracy');
legend('P-InceptionV3','P-Resnet50','P-Xception');
grid on;
tightfig();