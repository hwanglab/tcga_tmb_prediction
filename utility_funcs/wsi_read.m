%% this function is used to read RGB images with a certain resoluton

function RGB=wsi_read(slidePtr,objectivePower,downsampleFactors,width,height,magCoarse)

%% read the magCoarse image or less resolution for processing
mag=objectivePower./round(downsampleFactors);
xPos=1;yPos=1;

if any(mag==magCoarse)   %%there is magCoarse magnification
    levelToUse = find(mag == magCoarse,1);
    wid3=floor(width/downsampleFactors(levelToUse));
    hei3=floor(height/downsampleFactors(levelToUse));
    ARGB = openslide_read_region(slidePtr,xPos,yPos,wid3-1,hei3-1,levelToUse-1);
    RGB=ARGB(:,:,2:4);

else                     %% there is no magCoarse magnification in original image
    magToUseBelow=max(mag(mag<magCoarse));
    magToUseAbove=min(mag(mag>magCoarse));                  %% for reading image patches in a bit high resolution
    if isempty(magToUseBelow)                               %% there is not lower resolution than magCoarse
        levelToUse=find(mag==magToUseAbove,1);
        wid3=floor(width/downsampleFactors(levelToUse));
        hei3=floor(height/downsampleFactors(levelToUse));
        ARGB = openslide_read_region(slidePtr,xPos,yPos,wid3-1,hei3-1,levelToUse-1);
        RGB=ARGB(:,:,2:4);
        RGB=imresize(RGB,magCoarse/magToUseAbove);          %% reduce to magCoarse for preprocessing
    else
        levelToUse=find(mag == magToUseBelow,1);
        wid3=floor(width/downsampleFactors(levelToUse));
        hei3=floor(height/downsampleFactors(levelToUse));
        ARGB = openslide_read_region(slidePtr,xPos,yPos,wid3-1,hei3-1,levelToUse-1);
        RGB=ARGB(:,:,2:4);
        RGB=imresize(RGB,magCoarse/magToUseBelow);
    end
end