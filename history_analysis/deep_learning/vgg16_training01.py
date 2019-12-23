import os
from keras.models import Sequential
from keras.layers import Conv2D,MaxPool2D
from keras.layers import Activation,Dropout,Flatten,Dense
from keras.optimizers import Adam
from keras.preprocessing.image import ImageDataGenerator
from keras.callbacks import EarlyStopping,TensorBoard,ModelCheckpoint
import time

Train_data=['E:/Hongming/projects/tcga-bladder-mutationburden/deep_learning_data/data01/train/low/',
            'E:/Hongming/projects/tcga-bladder-mutationburden/deep_learning_data/data01/train/mid/',
            'E:/Hongming/projects/tcga-bladder-mutationburden/deep_learning_data/data01/train/high/']

Validation_data=['E:/Hongming/projects/tcga-bladder-mutationburden/deep_learning_data/data01/validation/low/',
                 'E:/Hongming/projects/tcga-bladder-mutationburden/deep_learning_data/data01/validation/mid/',
                 'E:/Hongming/projects/tcga-bladder-mutationburden/deep_learning_data/data01/validation/high/']

train_path='E:/Hongming/projects/tcga-bladder-mutationburden/deep_learning_data/data01/train/'
valid_path='E:/Hongming/projects/tcga-bladder-mutationburden/deep_learning_data/data01/validation/'

TENSORBOARD_PATH='E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/step8)deep_learning/assets/tensorboard/'
MODELS_PATH='E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/step8)deep_learning/assets/models/'


#Train_data=['E:/Hongming/deep_learning/keras_learning/data/train/cats/',
#            'E:/Hongming/deep_learning/keras_learning/data/train/dogs/',
#            'E:/Hongming/deep_learning/keras_learning/data/train/flowers/']

#Validation_data=['E:/Hongming/deep_learning/keras_learning/data/validation/cats/',
#            'E:/Hongming/deep_learning/keras_learning/data/validation/dogs/',
#            'E:/Hongming/deep_learning/keras_learning/data/validation/flowers/']

#train_path='E:/Hongming/deep_learning/keras_learning/data/train/'
#valid_path='E:/Hongming/deep_learning/keras_learning/data/validation/'

batch_size=16
sr=256
sc=256
def get_img_id():
    train_ids=[]
    for k in range(0,len(Train_data)):
        path_temp=Train_data[k]
        images=os.listdir(path_temp)
        for image_name in images:
            if '.png' or '.jpg' in image_name:
                train_ids.append(image_name)

    validation_ids=[]
    for j in range(0,len(Validation_data)):
        path_temp=Validation_data[j]
        images=os.listdir(path_temp)
        for image_name in images:
            if '.png' or '.jpg' in image_name:
                validation_ids.append(image_name)

    return train_ids, validation_ids

def get_Lecun_model():
    model=Sequential()
    model.add(Conv2D(32,(3,3),input_shape=(sr,sc,3))) # in input_shape: the batch dimension is not included
    model.add(Activation('relu'))
    model.add(MaxPool2D(pool_size=(2,2)))

    model.add(Conv2D(32,(3,3)))
    model.add(Activation('relu'))
    model.add(MaxPool2D(pool_size=(2,2)))

    model.add(Conv2D(64,(3,3)))
    model.add(Activation('relu'))
    model.add(MaxPool2D(pool_size=(2,2)))

    model.add(Flatten())
    model.add(Dense(64))
    model.add(Activation('relu'))
    model.add(Dropout(0.1))
    model.add(Dense(3))
    model.add(Activation('softmax'))

    model.compile(loss='categorical_crossentropy',optimizer=Adam(lr=1e-4),metrics=['accuracy'])
    return model

def train_classification():
    train_ids,validation_ids=get_img_id()
    model=get_Lecun_model()

    model_json = model.to_json()
    with open(MODELS_PATH + "model_first_try.json", "w") as json_file:
        json_file.write(model_json)
        print("save model to disk")

    train_datagen=ImageDataGenerator(
        rescale=1./255,
        shear_range=0.2,
        zoom_range=0.2,
        horizontal_flip=True)

    test_datagen=ImageDataGenerator(rescale=1./255)

    # Training new model
    ts = str(int(time.time()))
    model_name = 'VGG16_01'
    run_name = 'model={}-batch_size={}-ts={}'.format(model_name,batch_size,ts)
    tensorboard_loc = os.path.join(TENSORBOARD_PATH, run_name)
    checkpoint_loc = os.path.join(MODELS_PATH, 'model-{}-weights.h5'.format(ts))

    earlyStopping = EarlyStopping(monitor='val_loss',  # quantity to be monitored
                                  patience=5,
                                  # number of epochs with no improvement after which training will be stopped
                                  verbose=1,  # decides what to print
                                  min_delta=0.0001,  # minimum change to qualify as an improvement
                                  mode='min', )  # min mode: training will stop when the quantity monitored has stopped decreasing

    modelCheckpoint = ModelCheckpoint(checkpoint_loc,  # path to save the model
                                      monitor='val_loss',  # quantity to monitor
                                      save_best_only=True,  # if ture, the latest best model will not be overwritten
                                      mode='min',  # the decision to overwrite based on
                                      verbose=1,
                                      save_weights_only=True)

    tensorboard = TensorBoard(log_dir=tensorboard_loc, histogram_freq=0, write_graph=True, write_images=True)

    callbacks_list = [modelCheckpoint, earlyStopping, tensorboard]



    train_generator=train_datagen.flow_from_directory(
        directory=train_path,
        target_size=(sr,sc),
        batch_size=batch_size,
        class_mode='categorical')

    validation_generator=test_datagen.flow_from_directory(
        directory=valid_path,
        target_size=(sr,sc),
        batch_size=batch_size,
        class_mode='categorical'
    )

    model.fit_generator(
        train_generator,
        #callbacks=callbacks_list,
        steps_per_epoch=len(train_ids)//batch_size,
        epochs=30,
        verbose=1,
        validation_data=validation_generator,
        validation_steps=len(validation_ids)//batch_size)

    model.save_weights('frist_try.h5')
    # save the trained model

    

if __name__ == '__main__':
    train_classification()