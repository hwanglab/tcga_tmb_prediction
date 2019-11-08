function ff=xu_feature_bins_extractionIII(finput,fb)

%ff=zeros(1,2*size(finput,2)*fb);
ff_mean=[];
%ff_std=[];
for t=1:size(finput,2)
    ftemp=finput(:,t);
   [N,Y2] = histcounts(ftemp,fb);
    
    for k=1:length(Y2)-1
        lowv=Y2(k);
        highv=Y2(k+1);
        ind=find(ftemp>=lowv & ftemp<highv);
        %ff=[ff,mean(ftemp(ind)),std(ftemp(ind))]; 
        ff_mean=[ff_mean,mean(ftemp(ind))];
        %ff_std=[ff_std,std(ftemp(ind))];
    end
end

%ff=[ff_mean,ff_std];
ff=ff_mean;
