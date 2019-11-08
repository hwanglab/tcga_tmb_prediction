function ff=xu_feature_clustering_extractionII(finput,L)

%ff=zeros(1,2*size(finput,2)*fb);
ff=[];

for k=1:max(L(:))
    ind=find(L==k);
    ff_mean=[mean(finput(ind,:),1),length(ind)];
    ff=[ff,ff_mean];
end



