function [top_left_np,bottom_right_np,ind]=xu_patch_selection(np,top_left,bottom_right,slidePtr,levelforRead,magFine,magCoarse,magToUseAbove)

wb=0.85*255; % white background thrshold
if size(top_left,1)>np
    score=zeros(size(top_left,1),1);
    if nargin==7
        for tind=1:size(top_left,1)
            tlp=top_left(tind,:);
            brp=bottom_right(tind,:);
            tlp=(tlp-1).*(magFine/magCoarse)+1;
            brp=(brp-1).*(magFine/magCoarse)+1;
            ARGB = openslide_read_region(slidePtr,tlp(2),tlp(1),brp(2)-tlp(2),brp(1)-tlp(1),levelforRead-1);
            RGB=ARGB(:,:,2:4);
            
            [~,H,~]=normalizeStaining(RGB);
            Hr=H(:,:,1);
            bw=imbinarize(Hr,graythresh(Hr(Hr<wb)));
            bw=bwareaopen(~bw,20); % remove noisty pixels
            score(tind,1)=sum(bw(:));
        end
    else
        for tind=1:size(top_left,1)
            tlp=top_left(tind,:);
            brp=bottom_right(tind,:);
            tlp=(tlp-1).*(magToUseAbove/magCoarse)+1;
            brp=(brp-1).*(magToUseAbove/magCoarse)+1;
            ARGB = openslide_read_region(slidePtr,round(tlp(2)),round(tlp(1)),brp(2)-tlp(2),brp(1)-tlp(1),levelforRead-1);
            RGB=ARGB(:,:,2:4);
            RGB=imresize(RGB,magFine/magToUseAbove);
            
            [~,H,~]=normalizeStaining(RGB);
            Hr=H(:,:,1);
            bw=imbinarize(Hr,graythresh(Hr(Hr<wb)));
            bw=bwareaopen(~bw,20); % remove noisty pixels
            score(tind,1)=sum(bw(:));
        end
    end
    
    [temp,sortIndex]=sort(score,'descend');
    top_left_np=top_left(sortIndex(1:np),:);
    bottom_right_np=bottom_right(sortIndex(1:np),:);
    ind=sortIndex(1:np);
else
    top_left_np=top_left;
    bottom_right_np=bottom_right;
    ind=(1:size(top_left,1))';
end

