
import os
import numpy as np
from keras.preprocessing.image import ImageDataGenerator
from keras.models import Sequential
from keras.layers import Dropout, Flatten, Dense
from keras import applications
from keras.preprocessing import image
from keras.applications.resnet50 import preprocess_input, decode_predictions
from skimage.io import imshow
import scipy.io
import glob

# dimensions of image patches under 10x
img_width, img_height = 512, 512

#top_model_weights_path = 'bottleneck_fc_model.h5'
# train_data_dir = ['E:/Hongming/projects/tcga-bladder-mutationburden/nuclei_segmentation_imgs/mutation_analysis_imagesII/low/',
#                   'E:/Hongming/projects/tcga-bladder-mutationburden/nuclei_segmentation_imgs/mutation_analysis_imagesII/mid/',
#                   'E:/Hongming/projects/tcga-bladder-mutationburden/nuclei_segmentation_imgs/mutation_analysis_imagesII/high/']
#
# feat_output=['E:/Hongming/projects/tcga-bladder-mutationburden/featureoutput/7)transfer_learning/low/',
#              'E:/Hongming/projects/tcga-bladder-mutationburden/featureoutput/7)transfer_learning/mid/',
#              'E:/Hongming/projects/tcga-bladder-mutationburden/featureoutput/7)transfer_learning/high/']

train_data_dir =['E:/Hongming/projects/tcga-bladder-mutationburden/nuclei_segmentation_imgs/mutation_analysis_imagesIII/low/',
                 'E:/Hongming/projects/tcga-bladder-mutationburden/nuclei_segmentation_imgs/mutation_analysis_imagesIII/mid/',
                 'E:/Hongming/projects/tcga-bladder-mutationburden/nuclei_segmentation_imgs/mutation_analysis_imagesIII/high/']

feat_output=['E:/Hongming/projects/tcga-bladder-mutationburden/featureoutput/8)transfer_learning_all10x/low/',
             'E:/Hongming/projects/tcga-bladder-mutationburden/featureoutput/8)transfer_learning_all10x/mid/',
             'E:/Hongming/projects/tcga-bladder-mutationburden/featureoutput/8)transfer_learning_all10x/high/']

img_path=['E:/blca_mutationBurden/low/','E:/blca_mutationBurden/mid/','E:/blca_mutationBurden/high/']

def save_bottlebeck_features():
    #datagen = ImageDataGenerator(rescale=1. / 255)

    # build the resnet50 network
    model = applications.resnet50.ResNet50(weights='imagenet', include_top=False, pooling='avg')

    #img_train_ids=get_train_img_id()
    # load image setting the image size to 224 x 224
    for ind in range(0, len(img_path)):
        path_ind = img_path[ind]
        #path_split=path_ind.split("/")
        images = os.listdir(path_ind)
        for image_name in images:
            if '.svs' in image_name:
                patient_id=image_name[0:23]
                image_features = []
                image_names = []
                #patches=os.listdir(train_data_dir[ind]+patient_id+'*.png')
                patches=glob.glob(train_data_dir[ind]+patient_id+'*.png')

                for patch_name in patches:
                    patch_split=patch_name.split("\\")
                    img = image.load_img(patch_name, target_size=(512, 512))
                    # convert image to numpy array
                    x = image.img_to_array(img)
                    # the image is now in an array of shape (3, 224, 224)
                    # need to expand it to (1, 3, 224, 224) as it's expecting a list
                    x = np.expand_dims(x, axis=0)
                    x = preprocess_input(x)
                    # extract the features
                    features = model.predict(x)[0]
                    image_features.append(features)
                    image_names.append(patch_split[1])

                scipy.io.savemat(feat_output[ind]+patient_id+'_feat.mat', mdict={'image_features': image_features, 'image_names':image_names})
#
                # convert from Numpy to a list of values
                #features_arr = np.char.mod('%f', features)

    # generator = datagen.flow_from_directory(
    #     train_data_dir,
    #     target_size=(img_width, img_height), # dimension all images will be resized
    #     batch_size=batch_size,
    #     class_mode=None,
    #     shuffle=False)
    # bottleneck_features_train = model.predict_generator(
    #     generator, nb_train_samples // batch_size)
    # np.save(open('bottleneck_features_train.npy', 'wb'),
    #         bottleneck_features_train)
    #
    # generator = datagen.flow_from_directory(
    #     validation_data_dir,
    #     target_size=(img_width, img_height),
    #     batch_size=batch_size,
    #     class_mode=None,
    #     shuffle=False)
    # bottleneck_features_validation = model.predict_generator(
    #     generator, nb_validation_samples // batch_size)
    # np.save(open('bottleneck_features_validation.npy', 'wb'),
    #         bottleneck_features_validation)


# def train_top_model():
#     train_data = np.load(open('bottleneck_features_train.npy','rb'))
#     train_labels = np.array(
#         [0] * (nb_train_samples // 2) + [1] * (nb_train_samples // 2))
#
#     validation_data = np.load(open('bottleneck_features_validation.npy','rb'))
#     validation_labels = np.array(
#         [0] * (nb_validation_samples // 2) + [1] * (nb_validation_samples // 2))
#
#     model = Sequential()
#     model.add(Flatten(input_shape=train_data.shape[1:]))
#     model.add(Dense(256, activation='relu'))
#     model.add(Dropout(0.5))
#     model.add(Dense(1, activation='sigmoid'))
#
#     model.compile(optimizer='rmsprop',
#                   loss='binary_crossentropy', metrics=['accuracy'])
#
#     model.fit(train_data, train_labels,
#               epochs=epochs,
#               batch_size=batch_size,
#               validation_data=(validation_data, validation_labels))
#     model.save_weights(top_model_weights_path)
#
#     # save the trained model
#     model_json=model.to_json()
#     with open("model_bottleneckV01.json","w") as json_file:
#         json_file.write(model_json)
#         print("save model to disk")

save_bottlebeck_features()
#train_top_model()
