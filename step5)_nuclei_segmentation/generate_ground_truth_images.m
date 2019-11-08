%---creating nuclei ground truth mask from xml files---
% author: Hongming Xu, Cleveland Clinic


clc;
clearvars;
GTpath='E:\Hongming\resources\tcga_nuclei_segmentation\Annotations\';
%Imgpath='E:\Hongming\resources\tcga_nuclei_segmentation\Tissue images\';

Imgpath='E:\Hongming\resources\tcga_nuclei_segmentation\Hongming_Data\Normalized_TissueImgs\';

debugpath='E:\Hongming\resources\tcga_nuclei_segmentation\Hongming_Data\debug_path\';
outputGT='E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\ground_truth_images\';
list=dir(strcat(GTpath,'*.xml'));

img_type='.png';

SE=strel('disk',1);   % key parameter to erode nuclei boundary pixels
total=0;
total2=0;
for z=1:length(list)
    
    
    %image = imread(strcat('/Volumes/Surabhi/BTP/BTP_files/annotations/z/',tline));
    xml_file=strcat(GTpath,list(z).name);
    if exist(strcat(Imgpath,list(z).name(1:end-4),img_type))
        image=imread(strcat(Imgpath,list(z).name(1:end-4),img_type));
    else
        image=imread(strcat(Imgpath,list(z).name(1:end-4),'.tiff'));
    end
    %figure,imshow(image)
    [im_height, im_width, im_channel] = size(image);
    % return;
    count=0;
    xDoc = xmlread(xml_file);
    Regions=xDoc.getElementsByTagName('Region'); % get a list of all the region tags
    idd_total=[];
    for regioni = 0:Regions.getLength-1
        Region=Regions.item(regioni);  % for each region tag
        idd=Region.getAttribute('Id');
        idd_total=[idd_total,idd];
        verticies=Region.getElementsByTagName('Vertex'); %get a list of all the vertexes (which are in order)
        xy{regioni+1}=zeros(verticies.getLength-1,2); %allocate space for them
        for vertexi = 0:verticies.getLength-1 %iterate through all verticies
            x=str2double(verticies.item(vertexi).getAttribute('X')); %get the x value of that vertex
            y=str2double(verticies.item(vertexi).getAttribute('Y')); %get the y value of that vertex
            %hold on,plot(x,y,'y.')
           
            if(x<0)
                x=0;
            end
            if(y<0)
                y=0;
            end
            if(x>im_width)
                x=im_width;
            end
            if(y>im_height)
                y=im_height;
            end
            xy{regioni+1}(vertexi+1,:)=[x,y]; % finally save them into the array
            
        end
    end
    total2=total2+str2double(idd_total(end));
    %figure,hold all
    %set(gca,'YDir','reverse'); %invert y axis
    % for zz=1:length(xy)
    %     plot(xy{zz}(:,1),xy{zz}(:,2),'LineWidth',12)
    % end
    
    % svsinfo=imfinfo(svs_file);
    
    % s=1; %base level of maximum resolution
    % s2=2; % down sampling of 1:32
    %hratio=svsinfo(s2).Height/svsinfo(s).Height;  %determine ratio
    %wratio=svsinfo(s2).Width/svsinfo(s).Width;
   
    nrow=im_height;
    ncol=im_width;
    mask1=zeros(nrow,ncol); %pre-allocate a mask
    mask2=zeros(nrow,ncol); %pre-allocate a mask
    final=zeros(nrow,ncol,3); %pre-allocate a mask
    for zz=1:length(xy) %for each region
        count=count+1;
        smaller_x=xy{zz}(:,1); 
        smaller_y=xy{zz}(:,2);
        t = poly2mask(smaller_x,smaller_y,nrow,ncol);
        %t=imresize(t,0.5); % change to 20x resolution
        %t_d=imdilate(t,SE);
        %t_e=imerode(t,SE);   % the boundary on or within nuclei
        t_e=Iterative_erosion(t,SE);
        mask1=mask1+t-t_e;   % get the boundary regions
        mask2=mask2+t_e;     % get the nuclei regions
        
    end
    
    if count~=Regions.getLength
        disp('attention!!!!!!!!!!!');
    end
    clear xy xDoc Regions verticies t t_e Region
    total=total+count;
    
    mask1(mask1>0)=1;
    mask2(mask2>0)=1;
    final(:,:,1)=~(mask1 | mask2);  % background channel
    mask2=mask2-(mask2&mask1);      % make every pixel only belong one class
    final(:,:,2)=mask1;             % object channel
    final(:,:,3)=mask2;             % boundary channel
    % img(:,:,1)=final;
    % img(:,:,2)=final;
    % img(:,:,3)=final;
    imwrite(final,strcat(outputGT, strcat(list(z).name(1:end-4),'.png')));
    %img=imoverlay(image,mask1);
    %imwrite(img,strcat(debugpath,strcat(list(z).name(1:end-4),'.png')));
    clear mask1 mask2 final img
    
end
t=0;




