%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used for cal. the grid for given image tile size
% SELECT tiles only based on white foreground prior knowledge

%  Inputs:
%  -bwGT: the binary mask foreground (manually labeled)
%  -pp: the percent of the largest intrested value (the interest value is
%  the raito of nuclei pixels over the total pixels)
%  -tileSize: the image tile size

%  Output:
%  -top_left : each row is the top-left corner of a selected block; the
%  first column: row number; the second column: column number
%  -bottome_right: each row is the bottom-right corner of a selected block



% (c) Edited by Hongming Xu,
% Deptment of Quantitative Health Sciences,
% Cleveland Clinic, USA.  December 2017
% If you have any problem feel free to contact me.
% Please address questions or comments to: mxu@ualberta.ca

% Terms of use: You are free to copy,
% distribute, display, and use this work, under the following
% conditions. (1) You must give the original authors credit. (2) You may
% not use or redistribute this work for commercial purposes. (3) You may
% not alter, transform, or build upon this work. (4) For any reuse or
% distribution, you must make clear to others the license terms of this
% work. (5) Any of these conditions can be waived if you get permission
% from the authors.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# Usage: consider all the WSIs, but only selects all satisfied tiles
function [top_left,bottom_right]=xu_SelectImageTiles_VI(bwTissue,bwGT,pp,ppt,tileSize)

CC=bwconncomp(bwGT);
stats=regionprops(CC,'BoundingBox');
bb=cat(1,stats.BoundingBox);
ss=10; %% for safty not out of image, this parimeter is not affect performance
fun=@(x)sum(sum(x.data))/(tileSize(1)*tileSize(2));
top_left=[];
bottom_right=[];

for bi=1:size(bb,1)
    tbb=bb(bi,:);
    rs=round(tbb(2))+ss;
    cs=round(tbb(1))+ss;
    re=round(tbb(2))+round(tbb(4))-ss;
    ce=round(tbb(1))+round(tbb(3))-ss;
    Bgt=blockproc(bwGT(rs:re,cs:ce),tileSize,fun);        %% ground truth region
    Bt1=(Bgt>pp);
    BTis=blockproc(bwTissue(rs:re,cs:ce),tileSize,fun);   %% tissue binary mask
    Bt2=(BTis>ppt);
    B=Bt1&Bt2;
   
    yy=[rs,rs+tileSize(2):tileSize(2):re,re];           %% row direction
    xx=[cs,cs+tileSize(1):tileSize(1):ce,ce];           %% column direction
    [ry,cx]=find(B);
    top_left=[top_left;yy(ry)',xx(cx)'];              %% the first column: row; the second column: column; top-left point position
    bottom_right=[bottom_right;yy(ry+1)',xx(cx+1)'];  %% the first column: row; the second column: column; bottom-right point position
    
end

end