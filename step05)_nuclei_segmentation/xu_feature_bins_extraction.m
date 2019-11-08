function ff=xu_feature_bins_extraction(finput,fb)
temp=floor(100/fb);
bt=temp:temp:100-temp;

%ff=zeros(1,2*size(finput,2)*fb);
ff_mean=[];
ff_std=[];
for t=1:size(finput,2)
    ftemp=finput(:,t);
    Y=prctile(ftemp,bt);
    Y2=[min(ftemp),Y,max(ftemp)+0.001];
    for k=1:length(Y2)-1
        lowv=Y2(k);
        highv=Y2(k+1);
        ind=find(ftemp>=lowv & ftemp<highv);
        %ff=[ff,mean(ftemp(ind)),std(ftemp(ind))]; 
        ff_mean=[ff_mean,mean(ftemp(ind))];
       % ff_std=[ff_std,std(ftemp(ind))];
    end
end

ff=[ff_mean];
%ff=[ff_mean,ff_std];
