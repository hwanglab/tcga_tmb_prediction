'''
A simple example
step2: extract features

use pretrained classifiers as feature extractors

-- input: patches
-- output: feature vector saved into .mat files
author: Hongming Xu ccf,
email: mxu@ualberta.ca

reqirements: keras deep learning toolbox
'''
import os
import numpy as np
from keras.preprocessing.image import ImageDataGenerator
from keras.models import Sequential
from keras.layers import Dropout, Flatten, Dense, Input
from keras.preprocessing import image
from keras.applications import resnet50,nasnet,xception,inception_v3,inception_resnet_v2,densenet,vgg19
from keras.models import Model
from skimage.io import imshow
import scipy.io
import glob

# --- parameter settings ----#
# dimensions of image patches at 20x
img_width, img_height = 1024, 1024

wsi_dir='./'  # wsi patch
# selected tumor tiles by main_wsi_tiling_v01.m
data_dir='./step1_output/'

save_features=True # if True you will save features into .mat files
model_array=['resnet50','nasnet_large','xception','inceptionv3','inceptionresnetv2','densenet','vgg19']
# --- end parameter settings ----#

def save_bottlebeck_features(model_name):

    if model_name=='resnet50':
        model = resnet50.ResNet50(weights='imagenet', include_top=False, pooling='avg')

        # 2048 dimensional features
        # pooling: 1) None: output is 16x16x2048, 2) avg: 1x1x2048, 3) max: 1x1x2048
        #base_model=resnet50.ResNet50(input_shape=(img_height,img_width,3),weights='imagenet', include_top=False)
        #model = Model(inputs=base_model.input, outputs=base_model.get_layer('activation_25').output)
    elif model_name=='nasnet_large':
        model=nasnet.NASNetLarge(input_shape=(img_height,img_width,3),weights='imagenet', include_top=False, pooling='avg')
        #4032 dimensional features
    elif model_name=='xception':
        model=xception.Xception(input_shape=(img_height,img_width,3),weights='imagenet', include_top=False, pooling='avg')

        #2048 dimensional features
    elif model_name=='inceptionv3':
        model=inception_v3.InceptionV3(input_shape=(img_height,img_width,3),weights='imagenet', include_top=False, pooling='avg')

        #2048 dimensional features
    elif model_name=='inceptionresnetv2':
        model=inception_resnet_v2.InceptionResNetV2(input_shape=(img_height,img_width,3),weights='imagenet', include_top=False, pooling='avg')
        #1536 dimensional features
    elif model_name=='densenet':
        model=densenet.DenseNet201(input_shape=(img_height,img_width,3),weights='imagenet', include_top=False, pooling='avg')
        # 1920 dimensional features
    else:
        model=vgg19.VGG19(weights='imagenet', include_top=False, pooling='avg')
        # 512 dimensional features
        #base_model=vgg19.VGG19(input_shape=(img_height,img_width,3),weights='imagenet', include_top=False)
        #model=Model(inputs=base_model.input,outputs=base_model.get_layer('block4_pool').output)

    images = os.listdir(wsi_dir)
    for image_name in images:
        if '.svs' in image_name:
            patient_id=image_name[0:23]
            image_features = []
            image_names = []
            #patient_id='TCGA-2F-A9KT'
            #patches=os.listdir(train_data_dir[ind]+patient_id+'*.png')
            patches=glob.glob(data_dir+patient_id+'*.png')

            for patch_name in patches:
                patch_split=patch_name.split("\\")
                img = image.load_img(patch_name, target_size=(img_height,img_width))
                # convert image to numpy array
                x = image.img_to_array(img)

                # the image is now in an array of shape (224, 224, 3)
                # need to expand it to (1, 224, 224, 3) as it's expecting a list
                x = np.expand_dims(x, axis=0)
                #imshow(np.uint8(x[0,:,:,:]))

                if model_name=='resnet50':
                    x = resnet50.preprocess_input(x)
                elif model_name=='nasnet_large':
                    x = nasnet.preprocess_input(x)
                elif model_name == 'xception':
                    x = xception.preprocess_input(x)
                elif model_name=='inceptionv3':
                    x=inception_v3.preprocess_input(x)
                elif model_name == 'inceptionresnetv2':
                    x=inception_resnet_v2.preprocess_input(x)
                elif model_name=='densenet':
                    x=densenet.preprocess_input(x)
                else:
                    x=vgg19.preprocess_input(x)

                # extract the features
                features = model.predict(x)[0]
                #features=np.mean(features,axis=(0,1))

                image_features.append(features)
                image_names.append(patch_split[1])

            if save_features==True:
                scipy.io.savemat('./step2_output/'+patient_id+'_feat.mat', mdict={'image_features': image_features, 'image_names':image_names})


if __name__=='__main__':
    save_bottlebeck_features(model_array[2])

