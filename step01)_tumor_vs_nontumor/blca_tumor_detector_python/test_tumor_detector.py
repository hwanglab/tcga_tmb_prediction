'''
example for test tumor detector
input: wsi image
output: tumor prediction probability mask

author: Hongming Xu, CCF
email: mxu@ualberta.ca
'''

import os
import scipy.misc
import time
from keras.models import model_from_json
from wsi_tiling_prediction import wsi_tiling_prediction

wsi_path='E:/data/blca_mutationBurden/blca_wsi/'

model_weigth_path='./blca_tumor_detector/models/'

output_mask_path='./'

tile_size=512
magnification=20

def load_model_prediction():

    # step 1: load trained model
    # load json and create model
    json_file=open(model_weigth_path+'model_vgg16_like.json','r')
    loaded_model_json=json_file.read()
    json_file.close()
    loaded_model=model_from_json(loaded_model_json)

    #load weigths into new model
    loaded_model.load_weights(model_weigth_path+'model_vgg16_like_weights.h5')
    print("loaded model from disk")

    loaded_model.summary()
    # step 2: image tiling and prediction
    wsis = os.listdir(wsi_path)
    for img_name in wsis:
        if '.svs' in img_name:


            start_time=time.time()
            pmask=wsi_tiling_prediction(loaded_model,wsi_path+img_name, magnification, tile_size, MappingMag=2.5, Coverage=0.5, zscore=True)
            print("---{} minutes---".format((time.time()-start_time) / 60))
            scipy.misc.imsave(output_mask_path + img_name[0:23] + '.png', pmask)


if __name__=='__main__':
    load_model_prediction()