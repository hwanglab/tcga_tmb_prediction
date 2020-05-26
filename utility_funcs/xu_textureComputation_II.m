%%-----------usage-------------
% compute LBP texture features for image patches
% input? top_left-> top_left postions of image patches
%         bottom_right-> correspondingly bottom_right positions of image
%         patches
%         slidePtr-> pointer to image
%         levelforRead->level to read image
%         magFine-> rosolution to compute features
%         magCoarse-> rosolution for pre-processing
%         magToUseAbove->resultion higher than magFine
% output: feat: matrix with LBP featues, each row corresponds to a feature
% vector for each patch


function feat=xu_textureComputation_II(top_left,bottom_right,slidePtr,levelforRead,magFine,magCoarse,magToUseAbove)

%dd=1:2:5;                                       %% glcm distance
p=16;
mapping=getmapping(p,'riu2');
feat=zeros(size(top_left,1),(p+2)*4);

if nargin==6
    for tind=1:size(top_left,1)
        tlp=top_left(tind,:);
        brp=bottom_right(tind,:);
        tlp=(tlp-1).*(magFine/magCoarse)+1;
        brp=(brp-1).*(magFine/magCoarse)+1;
        ARGB = openslide_read_region(slidePtr,tlp(2),tlp(1),brp(2)-tlp(2),brp(1)-tlp(1),levelforRead-1);
        RGB=ARGB(:,:,2:4);
        
        feat(tind,:)=xu_LBP(RGB(:,:,1),mapping);
        
%         % b) second order glcm features
%         Hf=zeros(1,length(dd)*20);
%         for t=1:length(dd)
%             glcm = graycomatrix(RGB(:,:,1), 'offset', [0 dd(t);-dd(t) dd(t);-dd(t) 0;-dd(t) -dd(t)], 'Symmetric', true);
%             Har=xu_Haralick(glcm);   %<<20 Haralic features>>
%             Hf(1,20*(t-1)+1:20*t)=Har;
%         end
%         feat(tind,(p+2)*4+1:end)=Hf;
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
        
        feat(tind,:)=xu_LBP(RGB(:,:,1),mapping);
        
%         % b) second order glcm features
%         Hf=zeros(1,length(dd)*20);
%         for t=1:length(dd)
%             glcm = graycomatrix(RGB(:,:,1), 'offset', [0 dd(t);-dd(t) dd(t);-dd(t) 0;-dd(t) -dd(t)], 'Symmetric', true);
%             Har=xu_Haralick(glcm);   %<<20 Haralic features>>
%             Hf(1,20*(t-1)+1:20*t)=Har;
%         end
%         feat(tind,(p+2)*4+1:end)=Hf;
    end
end

