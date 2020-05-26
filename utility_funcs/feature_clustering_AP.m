%-------usage-----------
% perform AP clustering based on input features
% input: feat_tumor_cluster->feature matrix, each row corresponds to a
% feature vector
% output? idx->the clustering index for each sample
%          indf->the unique cluster index
% AP clustering reference:
% Frey, B. J., & Dueck, D. (2007). Clustering by passing messages between data points. science, 315(5814), 972-976.
% author: Hongming Xu, Ph.D., Cleveland Clinic
% feel free to use it, but should give original author credits


function [idx,indf]=feature_clustering_AP(feat_tumor_cluster)

debug=0;

N=size(feat_tumor_cluster,1);
x=feat_tumor_cluster;
M=N*N-N; s=zeros(M,3); % Make ALL N^2-N similarities
j=1;
for i=1:N
    for k=[1:i-1,i+1:N]
        %s(j,1)=i; s(j,2)=k; s(j,3)=-sum((x(i,:)-x(k,:)).^2);
        s(j,1)=i; s(j,2)=k; s(j,3)=-norm(x(i,:)-x(k,:));
        j=j+1;
    end
end
p=median(s(:,3)); % Set preference to median similarity
%p=min(s(:,3));
[idx,netsim,dpsim,expref]=apclusterSparse(s,p);
indf=unique(idx);

if debug==1
    fprintf('Number of clusters: %d\n',length(unique(idx)));
    fprintf('Fitness (net similarity): %f\n',netsim);
    %figure; % Make a figures showing the data and the clusters
    for i=unique(idx)'
        ii=find(idx==i); h=plot(x(ii,1),x(ii,2),'o'); hold on;
        col=rand(1,3); set(h,'Color',col,'MarkerFaceColor',col);
        xi1=x(i,1)*ones(size(ii)); xi2=x(i,2)*ones(size(ii));
        hold on,line([x(ii,1),xi1]',[x(ii,2),xi2]','Color',col);
    end
    axis equal tight;
end

