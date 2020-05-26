%----save selected representative patches by ap clustring--------%

function imgInfo=xu_save_selected_image_patches(top_left,bottom_right,slidePtr,levelforRead,magFine,magCoarse,imgname,imgpath,freq,HE,magToUseAbove)


imgInfo=cell(size(top_left,1),2);

if nargin==10
    for tind=1:size(top_left,1)
        tlp=top_left(tind,:);
        brp=bottom_right(tind,:);
        tlp=(tlp-1).*(magFine/magCoarse)+1;
        brp=(brp-1).*(magFine/magCoarse)+1;
        ARGB = openslide_read_region(slidePtr,tlp(2),tlp(1),brp(2)-tlp(2),brp(1)-tlp(1),levelforRead-1);
        RGB=ARGB(:,:,2:4);
        
        [RGB, ~, ~] = normalizeStaining_test(RGB,HE); % comment if P_E_CN
        
        if tind<10
            imgnn=strcat(imgname,num2str(0),num2str(0),num2str(0),num2str(tind),'.png');
            imwrite(RGB,strcat(imgpath,imgnn));
        elseif tind<100 && tind>=10
            imgnn=strcat(imgname,num2str(0),num2str(0),num2str(tind),'.png');
            imwrite(RGB,strcat(imgpath,imgnn));
        elseif tind<1000 && tind>=100
            imgnn=strcat(imgname,num2str(0),num2str(tind),'.png');
            imwrite(RGB,strcat(imgpath,imgnn));
        else
            imgnn=strcat(imgname,num2str(tind),'.png');
            imwrite(RGB,strcat(imgpath,imgnn));
        end
        
        imgInfo{tind,1}=imgnn;
        imgInfo{tind,2}=freq(tind);
        
    end
else
    for tind=1:size(top_left,1)
        tlp=top_left(tind,:);
        brp=bottom_right(tind,:);
        
        tlp=(tlp-1).*(magToUseAbove/magCoarse)+1;
        brp=(brp-1).*(magToUseAbove/magCoarse)+1;
        ARGB = openslide_read_region(slidePtr,tlp(2),tlp(1),brp(2)-tlp(2),brp(1)-tlp(1),levelforRead-1);
        RGB=ARGB(:,:,2:4);
        
        RGB=imresize(RGB,magFine/magToUseAbove);
        
        [RGB, ~, ~] = normalizeStaining_test(RGB,HE);  % P_E_CN
        
        if tind<10
            imgnn=strcat(imgname,num2str(0),num2str(0),num2str(0),num2str(tind),'.png');
            imwrite(RGB,strcat(imgpath,imgnn));
        elseif tind<100 && tind>=10
            imgnn=strcat(imgname,num2str(0),num2str(0),num2str(tind),'.png');
            imwrite(RGB,strcat(imgpath,imgnn));
        elseif tind<1000 && tind>=100
            imgnn=strcat(imgname,num2str(0),num2str(tind),'.png');
            imwrite(RGB,strcat(imgpath,imgnn));
        else
            imgnn=strcat(imgname,num2str(tind),'.png');
            imwrite(RGB,strcat(imgpath,imgnn));
        end
        
        imgInfo{tind,1}=imgnn;
        imgInfo{tind,2}=freq(tind);
        
    end
end

