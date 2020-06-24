'''
purpose: train tcga_blca_tumor_detector - to be used for detecting tumor regions in whole slide images
tumor detector architecture, see: https://www.biorxiv.org/content/10.1101/666503v1.abstract
author: Hongming Xu, 2019
email: mxu@ualberta.ca

To run the code, you needs:
python+keras+...
'''

import os
import time
from keras.models import Sequential
from keras.layers import Conv2D, MaxPool2D
from keras.layers import Activation, Dropout, Flatten, Dense
from keras.optimizers import Adam
from keras.preprocessing.image import ImageDataGenerator
import tensorflow as tf
import argparse
from keras.utils.training_utils import multi_gpu_model
from keras.preprocessing import image
import numpy as np
from keras import backend as K
from keras.callbacks import ModelCheckpoint,EarlyStopping,TensorBoard


# construct the argument parse and parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("-g", "--gpus", type=int, default=1,help="# of GPUs to use for training")
args = vars(ap.parse_args())

# grab the number of GPUs and store it in a conveience variable
G = args["gpus"]


Train_data=['E:/data_bladder_20x/train/non_tumor/',
            'E:/data_bladder_20x/train/tumor/']
Validation_data=['E:/data_bladder_20x/validation/non_tumor/',
                 'E:/data_bladder_20x/validation/tumor/']

train_path='E:/data_bladder_20x/train/'
valid_path='E:/data_bladder_20x/validation/'

TENSORBOARD_PATH='./tensorboard/'
MODELS_PATH='./models/'

batch_size = 128
rn=512
cn=512
kk=3
ps=3

def get_img_id():
    train_ids = []
    for k in range(0, len(Train_data)):
        path_temp = Train_data[k]
        images = os.listdir(path_temp)
        for image_name in images:
            if '.png' or '.jpg' in image_name:
                train_ids.append(image_name)

    validation_ids = []
    for j in range(0, len(Validation_data)):
        path_temp = Validation_data[j]
        images = os.listdir(path_temp)
        for image_name in images:
            if '.png' or '.jpg' in image_name:
                validation_ids.append(image_name)

    return train_ids, validation_ids


def get_VGG_Like_model():
    model = Sequential()
    model.add(Conv2D(8, (kk, kk), input_shape=(rn, cn, 3),padding='same',activation='relu'))  # in input_shape: the batch dimension is not included
    model.add(Conv2D(8, (kk, kk), padding='same', activation='relu'))
    model.add(MaxPool2D(pool_size=(ps, ps)))

    model.add(Conv2D(16, (kk, kk), padding='same', activation='relu'))
    model.add(Conv2D(16, (kk, kk), padding='same', activation='relu'))
    model.add(MaxPool2D(pool_size=(ps, ps)))

    model.add(Conv2D(32, (kk, kk),padding='same',activation='relu'))
    model.add(Conv2D(32, (kk, kk),padding='same',activation='relu'))
    model.add(Conv2D(32, (kk, kk),padding='same',activation='relu'))
    model.add(MaxPool2D(pool_size=(ps, ps)))

    model.add(Conv2D(64, (kk, kk),padding='same',activation='relu'))
    model.add(Conv2D(64, (kk, kk),padding='same',activation='relu'))
    model.add(Conv2D(64, (kk, kk),padding='same',activation='relu'))
    model.add(MaxPool2D(pool_size=(ps, ps)))

    model.add(Conv2D(64, (kk, kk),padding='same',activation='relu'))
    model.add(Conv2D(64, (kk, kk),padding='same',activation='relu'))
    model.add(Conv2D(64, (kk, kk),padding='same',activation='relu'))
    model.add(MaxPool2D(pool_size=(ps, ps)))

    model.add(Flatten())
    model.add(Dense(128,activation='relu'))
    model.add(Dropout(0.25))
    model.add(Dense(128,activation='relu'))
    model.add(Dropout(0.25))
    model.add(Dense(1))
    model.add(Activation('sigmoid'))

    model.compile(loss='binary_crossentropy', optimizer='rmsprop', metrics=['accuracy'])
    return model

def preprocessor(input):
    image_aug = image.img_to_array(input)
    for k in range(3):
        ac=np.random.uniform(0.9,1.1)
        bc=np.random.uniform(-10,10)
        image_aug[:,:,k]=ac*image_aug[:,:,k]+bc

    image_aug = np.clip(image_aug, 0, 255)
    image_aug=image.array_to_img(image_aug)
    return image_aug


def train_classification():
    train_ids, validation_ids = get_img_id()

    # check to see if we are compiling using just a single GPU
    if G <= 1:
        print("[INFO] training with 1 GPU...")
        model = get_VGG_Like_model()

    # otherwise, we are compiling using multiple GPUs
    else:
         print("[INFO] training with {} GPUs...".format(G))

         # we'll store a copy of the model on *every* GPU and then combine
         # the results from the gradient updates on the CPU
         with tf.device("/cpu:0"):
             # initialize the model
             model = get_VGG_Like_model()

         # make the model parallel
         model = multi_gpu_model(model, gpus=G)

    print(model.summary())


    train_datagen = ImageDataGenerator(
        # transformation first
        rotation_range=5,
        zoom_range=0.05,
        fill_mode='reflect',
        horizontal_flip=True,
        vertical_flip=True,
        channel_shift_range=20,
        brightness_range=[0.7, 1.3],

        # standardization second
        preprocessing_function=preprocessor,
        rescale=1 / 255,
        samplewise_center=True,
        samplewise_std_normalization=True)

    test_datagen = ImageDataGenerator(rescale=1. / 255,
                                      samplewise_center=True,
                                      samplewise_std_normalization=True)

    # Training new model
    ts = str(int(time.time()))
    model_name = 'VGG16_like'
    run_name = 'model={}-batch_size={}-ts={}'.format(model_name, batch_size, ts)
    tensorboard_loc = os.path.join(TENSORBOARD_PATH, run_name)
    checkpoint_loc = os.path.join(MODELS_PATH, 'model-{}-weights.h5'.format(ts))

    earlyStopping = EarlyStopping(monitor='val_loss',  # quantity to be monitored
                                  patience=5,          # number of epochs with no improvement after which training will be stopped
                                  verbose=1,           # decides what to print
                                  min_delta=0.0001,    # minimum change to qualify as an improvement
                                  mode='min', )        # min mode: training will stop when the quantity monitored has stopped decreasing

    modelCheckpoint = ModelCheckpoint(checkpoint_loc,       # path to save the model, save the model after every epoch
                                      monitor='val_loss',   # quantity to monitor
                                      save_best_only=True,  # if ture, the latest best model will not be overwritten
                                      mode='min',           # the decision to overwrite based on
                                      verbose=1,
                                      save_weights_only=True)

    tensorboard = TensorBoard(log_dir=tensorboard_loc, histogram_freq=0, write_graph=True, write_images=True)
    callbacks_list = [modelCheckpoint, earlyStopping, tensorboard]


    train_generator = train_datagen.flow_from_directory(
        directory=train_path,
        target_size=(rn, cn),
        batch_size=batch_size,
        class_mode='binary')

    validation_generator = test_datagen.flow_from_directory(
        directory=valid_path,
        target_size=(rn, cn),
        batch_size=batch_size,
        class_mode='binary'
    )

    model_json = model.to_json()
    with open(MODELS_PATH + "model_vgg16_like.json", "w") as json_file:
        json_file.write(model_json)
        print("save model to disk")


    model.fit_generator(
        train_generator,
        steps_per_epoch=len(train_ids) // batch_size,
        epochs=30,
        validation_data=validation_generator,
        validation_steps=len(validation_ids) // batch_size,
        callbacks=callbacks_list)

    #model.save_weights('frist_try.h5')
    # save the trained model
    #model_json = model.to_json()
    #with open("model_first_try.json", "w") as json_file:
    #    json_file.write(model_json)
    #    print("save model to disk")


if __name__ == '__main__':
    train_classification()
    print('training done!!')
