from __future__ import print_function

import os
from skimage.io import imsave, imread,imshow
from skimage import measure,filters, color, segmentation
import numpy as np
from keras.models import load_model
from keras import backend as K
from matplotlib import pyplot as plt
from training_nuclei_segmentation import ncce, w_categorical_crossentropy
from tqdm import tqdm
from scipy import ndimage

K.set_image_data_format('channels_last')  # TF dimension ordering in this code

MODELPATH='E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/step5)_nuclei_segmentation/Assets/models/'
TESTDATAPATH='E:/Hongming/projects/tcga-bladder-mutationburden/nuclei_segmentation_imgs/mutation_analysis_imagesII/high/'
pred_mask_dir='E:/Hongming/projects/tcga-bladder-mutationburden/nuclei_segmentation_imgs/mutation_analysis_masksII/high/'

#pred_edge_dir='E:/Hongming/projects/tcga-bladder-mutationburden/nuclei_segmentation_imgs/predict_edges/'

#TESTDATAPATH='E:/Hongming/projects/tcga-bladder-mutationburden/nuclei_segmentation_imgs/testing_images/'
#pred_mask_dir='E:/Hongming/projects/tcga-bladder-mutationburden/nuclei_segmentation_imgs/predict_masks/'


pr=pc=256
overlap=64

pt=0.4

def get_img_id():
    img_ids=[]
    images=os.listdir(TESTDATAPATH)
    for image_name in images:
        if '.png' in image_name:
            img_ids.append(image_name)

    return img_ids

def load_image_disk(img_id,folder=TESTDATAPATH):
    img=imread(os.path.join(folder,img_id),as_grey=False)
    return img

def generate_weight_map():
    pr2=pr+2  # to avaoid boundary 0
    pc2=pc+2
    bw_boundary=np.zeros((pr2,pc2),dtype=int)
    bw_boundary[1:-1,1:-1]=1
    bw_dis_boundary=ndimage.distance_transform_edt(bw_boundary)

    bw_center=np.ones((pr2,pc2),dtype=int)
    bw_center[int(pr2/2-1):int(pr2/2+1),int(pc2/2-1):int(pc2/2+1)]=0
    bw_dis_center=ndimage.distance_transform_edt(bw_center)

    ss=(pr2*pc2)/np.sum(np.divide(bw_dis_boundary,(bw_dis_boundary+bw_dis_center)))

    weight_map=ss*np.divide(bw_dis_boundary,(bw_dis_boundary+bw_dis_center))
    weight_map2=weight_map[1:-1,1:-1]
    return weight_map2



def image_segmentation(img_id,weight_map):
    img = load_image_disk(img_id)

    rr,cc=img.shape[:2]

    img_seg = np.zeros((rr, cc, 3), dtype=np.float32)
    img_coeff=np.zeros((rr,cc,3),dtype=np.float32)
    img_weight=np.dstack((weight_map,weight_map,weight_map))

    rs = 0
    re = rs + pr
    cs = 0
    ce = cs + pc

    while re <= rr:
        while ce <= cc:
            imgPatch = np.array(img[rs:re, cs:ce, :])
            imgPatch = np.array([imgPatch])
            imgs_mask_test = model.predict(imgPatch, verbose=1)
            img_seg[rs:re, cs:ce, :] = img_seg[rs:re, cs:ce, :] + np.multiply(np.squeeze(imgs_mask_test), img_weight)
            img_coeff[rs:re, cs:ce, :] = img_coeff[rs:re, cs:ce, :] + img_weight

            cs = cs + pc - overlap
            ce = cs + pc


        rs = rs + pr - overlap
        re = rs + pr
        cs = 0
        ce = cs + pc



    img_seg=np.divide(img_seg,img_coeff)

    # bwn=img_seg[:,:,1]>pt
    # edge_labels = measure.label(bwn, background=0)
    # img_edge = segmentation.mark_boundaries(img, edge_labels, color=(0, 1, 0))
    # imsave(os.path.join(pred_edge_dir, img_id + '_edge.png'), img_edge)

    img_seg=(img_seg*255).astype(np.uint8)
    imsave(os.path.join(pred_mask_dir,img_id+'_mask.png'),img_seg)


def test_img_segmentation():
    img_ids=get_img_id()

    weight_map = generate_weight_map()

    for img_id in tqdm(img_ids):
        image_segmentation(img_id,weight_map)




if __name__=='__main__':
    print('-' * 30)
    print('Loading saved model...')
    print('-' * 30)

    #model = load_model(os.path.join(MODELPATH, 'model-1535553973.h5'),
    #                   custom_objects={'dice_coef_loss': dice_coef_loss})
    model = load_model(os.path.join(MODELPATH, 'model-1536590144.h5'), custom_objects={'ncce': ncce})
    test_img_segmentation()