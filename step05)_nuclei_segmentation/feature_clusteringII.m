%%------%%
% input: featm: feature matrix
%        numc: the number of desired clusters
% output: find: index of selected features
function [idx,C01]=feature_clusteringII(featm,numc,methods)

featm00=bsxfun(@minus,featm,mean(featm));    % feature standardization
featm01=bsxfun(@rdivide,featm00,std(featm00));   % feature standardization


rng default; % for reproducibility

if strcmp(methods,'kmeans')
    opts=statset('Display','final');
    [idx,C]=kmeans(featm01,numc,'Replicates',5,'Options',opts);
elseif strcmp (methods,'gaussian')
    options = statset('MaxIter',300);
    GMModel = fitgmdist(featm01,numc,'RegularizationValue',0.01,'Options',options);
    idx=posterior(GMModel,featm01);
    C=GMModel.mu;
end

C00=bsxfun(@times,C,std(featm00));        % recover to original range
C01=bsxfun(@plus,C00,mean(featm));