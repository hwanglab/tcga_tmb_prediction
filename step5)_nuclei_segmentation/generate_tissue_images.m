Imgpath='E:\Hongming\resources\tcga_nuclei_segmentation\Tissue images\';
GTpath='E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\ground_truth_images\';

tissueimagepath='E:\Hongming\projects\tcga-bladder-mutationburden\nuclei_segmentation_imgs\tissue_images_normalized\';
debugpath='E:\Hongming\projects\tcga-bladder-mutationburden\debugoutput\nuclei\';
list=dir(strcat(Imgpath,'*f'));

for z=1:length(list)
    img_file=strcat(Imgpath,list(z).name);
    img=imread(img_file);
    %img=imresize(img,0.5);
    imwrite(img,strcat(tissueimagepath,list(z).name(1:23),'.png'));
    
    gt=imread(strcat(GTpath,list(z).name(1:23),'.png'));
    img2=imoverlay(img,gt(:,:,2));
    imwrite(img2,strcat(debugpath,list(z).name(1:23),'.png'));
end