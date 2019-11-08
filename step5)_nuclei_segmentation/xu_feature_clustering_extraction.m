function ff=xu_feature_clustering_extraction(finput,L)

%ff=zeros(1,2*size(finput,2)*fb);
ff_mean=[];
ff_std=[];

for k=1:max(L(:))
    ind=find(L==k);
    ff_mean=[ff_mean,mean(finput(ind,:),1)];
    ff_std=[ff_std,std(finput(ind,:),[],1)];
end


ff=[ff_mean,ff_std];
