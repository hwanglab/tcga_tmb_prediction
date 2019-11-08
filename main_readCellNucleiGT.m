addpath(genpath('E:\matlab_repository\toolboxes\openslide-matlab-master\'));
addpath(genpath('E:\matlab_repository\misc\'));
clear all;
xmlfolder='C:\Users\xuh3\Desktop\test\';
filePattern=sprintf('%s/*.xml',xmlfolder);
baseFilename=dir(filePattern);



magCoarse=40;
magMid=5.0;
magFine=20;

imagePath={'C:\Users\xuh3\Desktop\test\'};
for ip=1:length(imagePath)
    
    imgPath=imagePath{ip};
    
    imgs=dir(fullfile(imgPath,'*.svs'));
    
    
    for im=1:numel(imgs)
        file1=fullfile(imgPath,imgs(im).name);
        fprintf('filename=%s\n%d',file1,im);
        slidePtr=openslide_open(file1);
        [mppX,mppY,width,height,numberOfLevels,...
            downsampleFactors,objectivePower]=openslide_get_slide_properties(slidePtr);
        mag=objectivePower./round(downsampleFactors);
        
        
        %1) read magCoarse image
        RGB=wsi_read(slidePtr,objectivePower,downsampleFactors,width,height,magCoarse);
        
        %2) read xml lables
        tt=objectivePower/magCoarse;
        xmlname=strcat(imgs(im).name(1:end-3),'xml');
        fullFileName=fullfile(xmlfolder,xmlname);
        s=xml2struct(fullFileName);
        
        rr=s.Annotations.Annotation.Regions.Region;
        bwGT=false(size(RGB(:,:,1)));
        figure,imshow(RGB)
        if isstruct(rr)
            bw=false(size(RGB(:,:,1)));
            v=rr.Vertices;
            m=cell2mat(v.Vertex);
            temp=[m.Attributes];
            X={temp.X};
            X2=cell2mat(X');
            X3=str2num(X2);
            
            Y={temp.Y};
            Y2=cell2mat(Y');
            Y3=str2num(Y2);
            
            
            X3=round(X3/tt);
            Y3=round(Y3/tt);
            
            ind=sub2ind(size(bw),Y3,X3);
            bw(ind)=1;
            
            bw2=roipoly(bw,X3,Y3);
            bwGT=bwGT|bw2;
        else
            for n=1:length(rr)
                bw=false(size(RGB(:,:,1)));
                v=rr{1,n}.Vertices;
                x=str2double(v.Vertex.Attributes.X); %get the x value of that vertex
                y=str2double(v.Vertex.Attributes.Y); %get the y value of that vertex
               
                
                X3=round(x/tt);
                Y3=round(y/tt);
                
                ind=sub2ind(size(bw),Y3,X3);
                bw(ind)=1;
                
               
                bwGT=bwGT|bw;
                hold on,plot(X3,Y3,'g+');
                
            end
        end
        figure,imshow(bwGT);
        
    end
end

