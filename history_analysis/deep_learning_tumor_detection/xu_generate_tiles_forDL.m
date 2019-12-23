function xu_generate_tiles_forDL(top_left,bottom_right,slidePtr,levelforRead,magFine,magCoarse,imgname,imgpath,HE,magToUseAbove)

%imgoutput='E:\Hongming\projects\tcga-bladder-mutationburden\debugoutput\4)\';
% Wi=para.Wi;
% Hi=para.Hi;
% Hiv=para.Hiv;
% nstains=para.nstains;
% lambda=para.lambda;




if nargin==9
    for tind=1:size(top_left,1)
        tlp=top_left(tind,:);
        brp=bottom_right(tind,:);
        
        tlp=(tlp-1).*(magFine/magCoarse)+1;
        brp=(brp-1).*(magFine/magCoarse)+1;
        ARGB = openslide_read_region(slidePtr,tlp(2),tlp(1),brp(2)-tlp(2),brp(1)-tlp(1),levelforRead-1);
        RGB=ARGB(:,:,2:4);
        
        
        %[Wis, His,Hivs]=stainsep(RGB,nstains,lambda);
        %[RGB]=SCN(RGB,Hi,Wi,His);
        [RGB, ~, ~] = normalizeStaining_test(RGB,HE);
        
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
        
        
        %[Wis, His, Hivs]=stainsep(RGB,nstains,lambda);
        %[RGB]=SCN(RGB,Hi,Wi,His);
        
        [RGB, ~, ~] = normalizeStaining_test(RGB,HE);
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
        
    end
end

