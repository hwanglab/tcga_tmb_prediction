function ff=xu_feature_bins_extractionII(finput)


ff_mean=mean(finput,1);
ff_std=std(finput,[],1);
ff_skew=skewness(finput,[],1);
ff_kurtosis=kurtosis(finput,[],1);

ff=[ff_mean,ff_std,ff_skew,ff_kurtosis];
