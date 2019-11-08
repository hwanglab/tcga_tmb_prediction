import os
import numpy as np
from skimage.io import imsave, imread,imshow
from tqdm import tqdm
from PIL import Image

training_TissueImgs='E:/Hongming/projects/tcga-bladder-mutationburden/nuclei_segmentation_imgs/tissue_images_normalized/'
training_GTMasks='E:/Hongming/projects/tcga-bladder-mutationburden/nuclei_segmentation_imgs/ground_truth_images/'
Training_Input='E:/Hongming/projects/tcga-bladder-mutationburden/nuclei_segmentation_imgs/training_tissue_patches/'
Traing_GroundTruth='E:/Hongming/projects/tcga-bladder-mutationburden/nuclei_segmentation_imgs/training_gt_patches/'
pr=pc=256
overlap=70

def get_img_id():
    img_ids=[]
    images=os.listdir(training_TissueImgs)
    for image_name in images:
        #if '.tif' or '.tiff' in image_name:
        if '.png' in image_name:
            img_ids.append(image_name)

    mask_ids=[]
    masks=os.listdir(training_GTMasks)
    for mask_name in masks:
        if '.png' in mask_name:
            mask_ids.append(mask_name)

    return img_ids,mask_ids

def load_image_disk(img_id,folder=training_TissueImgs):
    img=imread(os.path.join(folder,img_id),as_grey=False)
    return img

def load_mask_disk(mask_id,folder=training_GTMasks):
    mask=imread(os.path.join(folder,mask_id),as_grey=False)
    return mask

def image_division(img_id):
    img = load_image_disk(img_id)
    if img_id.split('.')[1]=='png':
        mask=load_mask_disk(img_id[:-4]+'.png')
    elif img_id.split('.')[1]=='tiff':
        mask = load_mask_disk(img_id[:-5] + '.png')
    else:
        print('impossible!!!!')

    rr,cc=img.shape[:2]
    #yy=np.arange(0,rr,(pr-overlap))
    #yy=np.append(yy,rr)
    #xx=np.arange(0,cc,(pc-overlap))
    #xx=np.append(xx,cc)

    rs=0
    re=rs+pr
    cs=0
    ce=cs+pc
    r=c=0
    while re<=rr:
        while ce<=cc:
            imgPatch = np.array(img[rs:re, cs:ce, :])
        # im = Image.fromarray(imgPatch)
            imsave(Training_Input + img_id[:-4] + '_' + str(r) + '_' + str(c) + '.png', imgPatch)
        # imsave(Training_Input+img_id[:-4]+'_'+str(r)+'_'+str(c)+'.jpg',imgPatch,quality=100)

            maskPatch = np.array(mask[rs:re, cs:ce, :])
            imsave(Traing_GroundTruth + img_id[:-4] + '_' + str(r) + '_' + str(c) + '.png', maskPatch)

            cs = cs + pc - overlap
            ce = cs + pc
            c += 1

        rs=rs+pr-overlap
        re=rs+pr
        r+=1
        cs = 0
        ce = cs + pc


    # for r in range(len(yy)-1):
    #     rs=r*(pr-overlap)
    #     re=rs+pr
    #
    #     for c in range(len(xx)-1):
    #         cs=c*(pc-overlap)
    #         ce=cs+pc
    #
    #         imgPatch=np.array(img[rs:re,cs:ce,:])
    #         #im = Image.fromarray(imgPatch)
    #         imsave(Training_Input+img_id[:-4]+'_'+str(r)+'_'+str(c)+'.png',imgPatch)
    #         #imsave(Training_Input+img_id[:-4]+'_'+str(r)+'_'+str(c)+'.jpg',imgPatch,quality=100)
    #
    #         maskPatch=np.array(mask[rs:re,cs:ce,:])
    #         imsave(Traing_GroundTruth + img_id[:-4] + '_' + str(r) + '_' + str(c) + '.png', maskPatch)

def train_patch_generation():
    img_ids,mask_ids=get_img_id()

    for img_id in tqdm(img_ids):
        image_division(img_id)




if __name__=='__main__':
    train_patch_generation()