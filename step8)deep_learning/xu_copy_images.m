function xu_copy_images(source,destination,images)
for i=1:length(images)
    img_name=images(i).name;
    copyfile(strcat(source,img_name),strcat(destination,img_name));
end