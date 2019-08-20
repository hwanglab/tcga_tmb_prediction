
% for usage refere to: https://www.mathworks.com/help/stats/savecompactmodel.html
%saveCompactModel(trainedModel_SVM_Cubic.ClassificationSVM,'SVM_Cubic');

%saveCompactModel(trainedModel_SVM_FineGaussian.ClassificationSVM,'SVM_FineGaussian');

save('SVM_cubic_model.mat','trainedModel_SVM_cubic')
save('SVM_FineGaussian_model.mat','trainedModel_SVM_fine_Gaussian')