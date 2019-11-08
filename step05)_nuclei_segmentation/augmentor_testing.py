import Augmentor as Aug

#for k in range(0,5):
p=Aug.Pipeline("E:/Hongming/projects/tcga-bladder-mutationburden/nuclei_segmentation_imgs/augmentor_testing/")
    #p.ground_truth("E:/Hongming/resources/tcga_nuclei_segmentation/Hongming_Data/Training_GTMasks/")

#p.rotate(probability=0.8, max_left_rotation=10, max_right_rotation=10)         # random rotation

#p.flip_left_right(probability=0.8)                                           # random flip
#p.flip_top_bottom(probability=0.8)

p.random_distortion(probability=0.5, grid_width=8, grid_height=8, magnitude=5) # randome elastic transformation

#p.zoom(probability=0.3, min_factor=1.1, max_factor=1.6)


    #p.rotate(probability=0.5, max_left_rotation=10, max_right_rotation=10)
#p.zoom(probability=0.8, min_factor=1.1, max_factor=1.5)
    #p.skew(probability=0.8)
p.random_contrast(probability=0.8,min_factor=0.9,max_factor=1.4)
#p.sample(30)
p.process()   # process each image in pipeline once
