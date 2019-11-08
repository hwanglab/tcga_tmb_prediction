
import os
import numpy as np
import time
import pickle
from keras.models import Sequential
from keras.layers import Dense, Activation, Dropout, Flatten
from keras.layers import Conv2D
from keras.layers import MaxPooling2D
from skimage.io import imsave, imread,imshow
from keras.utils import to_categorical
from keras.callbacks import EarlyStopping, ModelCheckpoint, TensorBoard
from keras.optimizers import Adam, rmsprop
from keras import regularizers


train_path=['E:/Hongming/projects/tcga-bladder-mutationburden/tumor_patches_5x/low/',
           'E:/Hongming/projects/tcga-bladder-mutationburden/tumor_patches_5x/mid/',
           'E:/Hongming/projects/tcga-bladder-mutationburden/tumor_patches_5x/high/']

valid_path=['E:/Hongming/projects/tcga-bladder-mutationburden/tumor_patches_5x/validation/low/',
            'E:/Hongming/projects/tcga-bladder-mutationburden/tumor_patches_5x/validation/mid/',
            'E:/Hongming/projects/tcga-bladder-mutationburden/tumor_patches_5x/validation/high/']

TENSORBOARD_PATH='E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/step4)_wsi_classification/Assets/vggtensorboards/'
MODELS_PATH='E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/step4)_wsi_classification/Assets/vggmodels/'

input_shape = (256, 256, 3)
num_classes=3
BATCH_SIZE=32
ff_channel=16
num_epochs = 80

def get_train_img_id():
    img_train_ids = []
    img_train_labels={}
    for ind in range(0, len(train_path)):
        path_ind = train_path[ind]
        path_split=path_ind.split("/")
        images = os.listdir(path_ind)
        for image_name in images:
            if '.jpg' in image_name:
                img_train_ids.append(image_name)
                img_train_labels[image_name]=path_split[-2]

    return img_train_ids,img_train_labels

def get_valid_img_id():
    img_valid_ids = []
    img_valid_labels={}
    for ind in range(0, len(valid_path)):
        path_ind = valid_path[ind]
        path_split=path_ind.split("/")
        images = os.listdir(path_ind)
        for image_name in images:
            if '.jpg' in image_name:
                img_valid_ids.append(image_name)
                img_valid_labels[image_name]=path_split[-2]

    return img_valid_ids,img_valid_labels

def generate_training_batch(data, batch_size,img_labels):
    while True:
        X_batch = []
        Y_batch = []
        batch_ids = np.random.choice(data,size=batch_size,replace=False)
        for idx, img_id in enumerate(batch_ids):
            loc_label=img_labels[img_id]
            if loc_label=='low':
                temp_path=train_path[0]
                y=0
            elif loc_label=='mid':
                temp_path=train_path[1]
                y=1
            elif loc_label=='high':
                temp_path=train_path[2]
                y=2
            else:
                temp_path=0
                y=-1
                print("impossible~~~~~")

            x = imread(os.path.join(temp_path, img_id), as_grey=False)

            X_batch.append(x)
            Y_batch.append(y)

        X = np.asarray(X_batch, dtype=np.float32)
        Y = np.asarray(Y_batch, dtype=np.float32)
        Y2= to_categorical(Y, num_classes)
        #return X,Y2
        yield X, Y2

def generate_validation_batch(data, batch_size,img_labels):
    while True:
        X_batch = []
        Y_batch = []
        batch_ids = np.random.choice(data,size=batch_size, replace=False)
        for idx, img_id in enumerate(batch_ids):
            loc_label = img_labels[img_id]
            if loc_label == 'low':
                temp_path = valid_path[0]
                y = 0
            elif loc_label == 'mid':
                temp_path = valid_path[1]
                y = 1
            elif loc_label == 'high':
                temp_path = valid_path[2]
                y = 2
            else:
                temp_path = 0
                y = -1
                print("impossible~~~~~")

            x = imread(os.path.join(temp_path, img_id), as_grey=False)

            X_batch.append(x)
            Y_batch.append(y)

        X = np.asarray(X_batch, dtype=np.float32)
        Y = np.asarray(Y_batch, dtype=np.float32)
        Y2 = to_categorical(Y, num_classes)
        yield X, Y2
        #return X, Y2

def get_vgg():
    reg = regularizers.l1_l2(l1=0.01, l2=0.01)
    model = Sequential([
        Conv2D(ff_channel*1, (3, 3), input_shape=input_shape, padding='same',activation='selu'),
        Conv2D(ff_channel*1, (3, 3), activation='selu', padding='same'),
        MaxPooling2D(pool_size=(2, 2), strides=(2, 2)),

        #Dropout(0.25),
        Conv2D(ff_channel*2, (3, 3), activation='selu', padding='same'),
        Conv2D(ff_channel*2, (3, 3), activation='selu', padding='same',),
        MaxPooling2D(pool_size=(2, 2), strides=(2, 2)),

        #Dropout(0.25),
        Conv2D(ff_channel*4, (3, 3), activation='selu', padding='same',),
        Conv2D(ff_channel*4, (3, 3), activation='selu', padding='same',),
        Conv2D(ff_channel*4, (3, 3), activation='selu', padding='same',),
        MaxPooling2D(pool_size=(2, 2), strides=(2, 2)),

        #Dropout(0.25),
        Conv2D(ff_channel*8, (3, 3), activation='selu', padding='same',),
        Conv2D(ff_channel*8, (3, 3), activation='selu', padding='same',),
        Conv2D(ff_channel*8, (3, 3), activation='selu', padding='same',),
        MaxPooling2D(pool_size=(2, 2), strides=(2, 2)),

        #Dropout(0.25),
        Conv2D(ff_channel*8, (3, 3), activation='selu', padding='same',),
        Conv2D(ff_channel*4, (3, 3), activation='selu', padding='same',),
        Conv2D(ff_channel*2, (3, 3), activation='selu', padding='same',),
        MaxPooling2D(pool_size=(2, 2), strides=(2, 2)),

        #Dropout(0.25),
        Flatten(),
        #Dense(2048, activation='selu'),
        Dense(1024, activation='selu',W_regularizer=reg),
        Dense(num_classes, activation='softmax')
        ])

    model.compile(loss='categorical_crossentropy',
                  #optimizer=Adam(lr=1e-5),
                  optimizer='rmsprop',
                  metrics=['accuracy'])
    return model

if __name__ == '__main__':
    img_train_ids, img_train_labels=get_train_img_id()
    img_valid_ids, img_valid_labels = get_valid_img_id()

    model=get_vgg()
    print(model.summary())

    train_generator = generate_training_batch(img_train_ids, BATCH_SIZE,img_train_labels)
    valid_generator = generate_validation_batch(img_valid_ids, BATCH_SIZE,img_valid_labels)

    # Training new model
    ts = str(int(time.time()))
    model_name = 'vgg16'
    steps_per_epoch = int(len(img_train_ids ) / BATCH_SIZE)
    run_name = 'model={}-num_epoch={}-steps_per_epoch={}-ts={}'.format(model_name,num_epochs,steps_per_epoch,ts)
    tensorboard_loc = os.path.join(TENSORBOARD_PATH, run_name)
    checkpoint_loc = os.path.join(MODELS_PATH, 'model-{}-weights.h5'.format(ts))

    earlyStopping = EarlyStopping(monitor='val_loss',  # quantity to be monitored
                                  patience=10,
                                  # number of epochs with no improvement after which training will be stopped
                                  verbose=1,  # decides what to print
                                  min_delta=0.00001,  # minimum change to qualify as an improvement
                                  mode='min', )  # min mode: training will stop when the quantity monitored has stopped decreasing

    modelCheckpoint = ModelCheckpoint(checkpoint_loc,  # path to save the model
                                      monitor='val_loss',  # quantity to monitor
                                      save_best_only=True,  # if ture, the latest best model will not be overwritten
                                      mode='min',  # the decision to overwrite based on
                                      verbose=1,
                                      save_weights_only=True)

    tensorboard = TensorBoard(log_dir=tensorboard_loc, histogram_freq=0, write_graph=True, write_images=True)

    callbacks_list = [modelCheckpoint, earlyStopping, tensorboard]

    VALIDATION_STEPS = int(len(img_valid_ids) / BATCH_SIZE)

    print('Starting run {}'.format(run_name))
    history = model.fit_generator(
        train_generator,
        steps_per_epoch=steps_per_epoch,
        epochs=num_epochs,
        callbacks=callbacks_list,
        verbose=1,
        validation_data=valid_generator,
        validation_steps=VALIDATION_STEPS)

    model_path = os.path.join(MODELS_PATH, 'model-{}.h5'.format(ts))
    history_path = os.path.join(MODELS_PATH, 'model-{}.history'.format(ts))
    model.save(model_path)
    pickle.dump(history.history, open(history_path, "wb"))
    print('Saved model at {}'.format(model_path))
    print('Saved model history at {}'.format(history_path))