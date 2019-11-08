import os
import time
import numpy as np
from skimage.io import imsave, imread,imshow
from tqdm import tqdm
from keras.callbacks import EarlyStopping, ModelCheckpoint, TensorBoard
from keras.models import Model
from keras.layers import Input, concatenate, Conv2D, MaxPooling2D, Conv2DTranspose, Dropout
from keras.optimizers import Adam
from keras import backend as K
from sklearn import model_selection, preprocessing, metrics
import cv2
import matplotlib.pyplot as plt
import Augmentor as Aug
import PIL
import shutil
import pickle
import functools
from itertools import product


Training_Input='E:/Hongming/projects/tcga-bladder-mutationburden/nuclei_segmentation_imgs/training_tissue_patches/'
Training_GroundTruth='E:/Hongming/projects/tcga-bladder-mutationburden/nuclei_segmentation_imgs/training_gt_patches/'
MODELS_PATH = 'E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/step5)_nuclei_segmentation/Assets/models/'
TENSORBOARD_PATH = 'E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/step5)_nuclei_segmentation/Assets/tensorboard/'
temp_destination_Input='E:/Hongming/projects/tcga-bladder-mutationburden/nuclei_segmentation_imgs/temp_Input/'
temp_destination_GT='E:/Hongming/projects/tcga-bladder-mutationburden/nuclei_segmentation_imgs/temp_GT/'
temp_input_augmented='E:/Hongming/projects/tcga-bladder-mutationburden/nuclei_segmentation_imgs/temp_Input/output/'

HEIGHT = 256
WIDTH = 256
CHANNELS = 3
BATCH_SIZE = 8
num_epochs = 80
drop_pp=0.2
thr_pp=0.12

w_array = np.ones((3,3))


def get_img_id():
    img_ids=[]
    images=os.listdir(Training_Input)
    for image_name in images:
        if '.png' in image_name:
            img_ids.append(image_name)

    mask_ids=[]
    masks=os.listdir(Training_GroundTruth)
    for mask_name in masks:
        if '.png' in mask_name:
            mask_ids.append(mask_name)

    return img_ids,mask_ids

def get_img_id_xu(): # only load color images
    img_ids=[]
    images=os.listdir(temp_input_augmented)
    for image_name in images:
        if 'temp_Input_original_' in image_name:
            img_ids.append(image_name)
    return img_ids

def load_image_disk(img_id,folder=Training_Input):
    img=imread(os.path.join(folder,img_id),as_grey=False)
    return img

def load_mask_disk(mask_id,folder=Training_GroundTruth):
    mask=imread(os.path.join(folder,mask_id),as_grey=False)
    return mask

def dice_coef(y_true, y_pred):
    y_true_f0 = K.flatten(y_true)
    y_pred_f0 = K.flatten(y_pred)
    intersection0 = K.sum(y_true_f0 * K.log(y_pred_f0+1e-9))
    return intersection0


def dice_coef_loss(y_true, y_pred):
    return -dice_coef(y_true, y_pred)

# def w_categorical_crossentropy(y_true, y_pred, weights):
#     nb_cl = len(weights)
#     y_true=K.reshape(y_true,[-1,CHANNELS])
#     y_pred=K.reshape(y_pred,[-1,CHANNELS])
#     final_mask = K.zeros_like(y_pred[:, 0])
#     y_pred_max = K.max(y_pred, axis=1)
#     y_pred_max = K.expand_dims(y_pred_max,1)
#     y_pred_max_mat = K.equal(y_pred, y_pred_max)
#     for c_p, c_t in product(range(nb_cl), range(nb_cl)):
#         final_mask += (K.cast(weights[c_t, c_p],K.floatx()) * K.cast(y_pred_max_mat[:, c_p],K.floatx()) * K.cast(y_true[:, c_t],K.floatx()))
#     #return K.sum(K.categorical_crossentropy(y_pred, y_true) * final_mask)
#     return - K.sum(K.sum(y_true * K.log(y_pred+1e-9))*final_mask)

def w_categorical_crossentropy(y_true, y_pred, weights):
    nb_cl = len(weights)
    final_mask = K.zeros_like(y_pred[..., 0])
    y_pred_max = K.max(y_pred, axis=-1)
    y_pred_max = K.expand_dims(y_pred_max, axis=-1)
    y_pred_max_mat = K.equal(y_pred, y_pred_max)
    for c_p, c_t in product(range(nb_cl), range(nb_cl)):
        w = K.cast(weights[c_t, c_p], K.floatx())
        y_p = K.cast(y_pred_max_mat[..., c_p], K.floatx())
        y_t = K.cast(y_true[..., c_t], K.floatx())
        final_mask += w * y_p * y_t
    return K.categorical_crossentropy(y_true, y_pred) * final_mask

def ncce(y_true, y_pred):
    val=w_categorical_crossentropy(y_true, y_pred, weights=w_array)
    return val

#ncce = functools.partial(w_categorical_crossentropy, weights=w_array)
#ncce.__name__ ='w_categorical_crossentropy'

def get_unet():
    inputs = Input(shape=(HEIGHT, WIDTH, CHANNELS))  # change into 3D rgb image
    conv1 = Conv2D(32, (3, 3), activation='selu', padding='same')(inputs)
    #drop1=Dropout(drop_pp)(conv1)
    conv1 = Conv2D(32, (3, 3), activation='selu', padding='same')(conv1)
    #drop1=Dropout(drop_pp)(conv1)
    pool1 = MaxPooling2D(pool_size=(2, 2))(conv1)

    conv2 = Conv2D(64, (3, 3), activation='selu', padding='same')(pool1)
    #drop2=Dropout(drop_pp)(conv2)
    conv2 = Conv2D(64, (3, 3), activation='selu', padding='same')(conv2)
    #drop2=Dropout(drop_pp)(conv2)
    pool2 = MaxPooling2D(pool_size=(2, 2))(conv2)

    conv3 = Conv2D(128, (3, 3), activation='selu', padding='same')(pool2)
    #drop3=Dropout(drop_pp)(conv3)
    conv3 = Conv2D(128, (3, 3), activation='selu', padding='same')(conv3)
    #drop3 = Dropout(drop_pp)(conv3)
    pool3 = MaxPooling2D(pool_size=(2, 2))(conv3)

    conv4 = Conv2D(256, (3, 3), activation='selu', padding='same')(pool3)
    #drop4 = Dropout(drop_pp)(conv4)
    conv4 = Conv2D(256, (3, 3), activation='selu', padding='same')(conv4)
    #drop4 = Dropout(drop_pp)(conv4)
    pool4 = MaxPooling2D(pool_size=(2, 2))(conv4)

    conv5 = Conv2D(512, (3, 3), activation='selu', padding='same')(pool4)
    #drop5 = Dropout(drop_pp)(conv5)
    conv5 = Conv2D(512, (3, 3), activation='selu', padding='same')(conv5)


    up6 = concatenate([Conv2DTranspose(256, (2, 2), strides=(2, 2), padding='same')(conv5), conv4], axis=3)
    conv6 = Conv2D(256, (3, 3), activation='selu', padding='same')(up6)
    #drop6 = Dropout(drop_pp)(conv6)
    conv6 = Conv2D(256, (3, 3), activation='selu', padding='same')(conv6)

    up7 = concatenate([Conv2DTranspose(128, (2, 2), strides=(2, 2), padding='same')(conv6), conv3], axis=3)
    conv7 = Conv2D(128, (3, 3), activation='selu', padding='same')(up7)
    #drop7 = Dropout(drop_pp)(conv7)
    conv7 = Conv2D(128, (3, 3), activation='selu', padding='same')(conv7)

    up8 = concatenate([Conv2DTranspose(64, (2, 2), strides=(2, 2), padding='same')(conv7), conv2], axis=3)
    conv8 = Conv2D(64, (3, 3), activation='selu', padding='same')(up8)
    #drop8 = Dropout(drop_pp)(conv8)
    conv8 = Conv2D(64, (3, 3), activation='selu', padding='same')(conv8)

    up9 = concatenate([Conv2DTranspose(32, (2, 2), strides=(2, 2), padding='same')(conv8), conv1], axis=3)
    conv9 = Conv2D(32, (3, 3), activation='selu', padding='same')(up9)
    #drop9 = Dropout(drop_pp)(conv9)
    conv9 = Conv2D(32, (3, 3), activation='selu', padding='same')(conv9)

    conv10 = Conv2D(3, (1, 1), activation='softmax')(conv9)

    model = Model(inputs=[inputs], outputs=[conv10])

    #model.compile(optimizer=Adam(lr=1e-5), loss=dice_coef_loss, metrics=['accuracy'])
    #model.compile(optimizer=Adam(lr=1e-5), loss='categorical_crossentropy', metrics=['accuracy'])
    model.compile(optimizer=Adam(lr=1e-5), loss=ncce, metrics=['accuracy'])

    return model

def get_model_memory_usage(batch_size, model):

    shapes_mem_count = 0
    for l in model.layers:
        single_layer_mem = 1
        for s in l.output_shape:
            if s is None:
                continue
            single_layer_mem *= s
        shapes_mem_count += single_layer_mem

    trainable_count = int(np.sum([K.count_params(p) for p in set(model.trainable_weights)]))
    non_trainable_count = int(np.sum([K.count_params(p) for p in set(model.non_trainable_weights)]))

    total_memory = 4 * batch_size * (shapes_mem_count + trainable_count + non_trainable_count)
    gbytes = round(total_memory / (1024 ** 3), 3)
    mbytes = round(total_memory / (1024 ** 2), 3)

    print('trainable_count', trainable_count, 'non_trainable_count', non_trainable_count, 'gbytes', gbytes, 'mbytes',
          mbytes)

def randomShiftScaleRotate(image, mask,
                           shift_limit=(-0.0625, 0.0625),
                           scale_limit=(-0.1, 0.1),
                           rotate_limit=(-45, 45), aspect_limit=(0, 0),
                           borderMode=cv2.BORDER_REFLECT_101, u=0.5):
    if np.random.random() < u:
        height, width, channel = image.shape

        angle = np.random.uniform(rotate_limit[0], rotate_limit[1])  # degree
        scale = np.random.uniform(1 + scale_limit[0], 1 + scale_limit[1])
        aspect = np.random.uniform(1 + aspect_limit[0], 1 + aspect_limit[1])
        sx = scale * aspect / (aspect ** 0.5)
        sy = scale / (aspect ** 0.5)
        dx = round(np.random.uniform(shift_limit[0], shift_limit[1]) * width)
        dy = round(np.random.uniform(shift_limit[0], shift_limit[1]) * height)

        cc = np.math.cos(angle / 180 * np.math.pi) * sx
        ss = np.math.sin(angle / 180 * np.math.pi) * sy
        rotate_matrix = np.array([[cc, -ss], [ss, cc]])

        box0 = np.array([[0, 0], [width, 0], [width, height], [0, height], ])
        box1 = box0 - np.array([width / 2, height / 2])
        box1 = np.dot(box1, rotate_matrix.T) + np.array([width / 2 + dx, height / 2 + dy])

        box0 = box0.astype(np.float32)
        box1 = box1.astype(np.float32)
        mat = cv2.getPerspectiveTransform(box0, box1)
        image = cv2.warpPerspective(image, mat, (width, height), flags=cv2.INTER_LINEAR, borderMode=borderMode,
                                    borderValue=(0, 0, 0,))
        mask2 = cv2.warpPerspective(mask, mat, (width, height), flags=cv2.INTER_LINEAR, borderMode=borderMode,
                                   borderValue=(0, 0, 0,))

        # plt.figure(num='debug',figsize=(8,8))
        # plt.subplot(2,2,1)
        # plt.imshow(image)
        #
        # plt.subplot(2,2,2)
        # imshow(image2)
        #
        # plt.subplot(2,2,3)
        # plt.imshow(mask)
        #
        # plt.subplot(2,2,4)
        # plt.imshow(mask2)
        #
        # plt.show()


        ind_max=np.argmax(mask2,axis=2)
        mask2=np.zeros((HEIGHT,WIDTH,CHANNELS),np.uint8)
        x = np.linspace(0, HEIGHT-1, HEIGHT, dtype=np.int)
        y = np.linspace(0, WIDTH-1, WIDTH, dtype=np.int)
        xv, yv = np.meshgrid(x, y)
        mask2[yv,xv,ind_max]=255
        if len(mask.shape) == 2:
            mask = np.expand_dims(mask, axis=2)

    return image, mask

def generate_training_batch(data, batch_size,train_imgs,train_masks):
    while True:
        X_batch = []
        Y_batch = []
        batch_ids = np.random.choice(data,
                                     size=batch_size,
                                     replace=False)
        for idx, img_id in enumerate(batch_ids):
            #x = get_image(img_id)
            x=train_imgs[img_id]
            mask_id=img_id[:-3]+'png'
            #y = get_mask(mask_id)
            y=train_masks[mask_id]

            x, y = randomShiftScaleRotate(x, y,
                                          shift_limit=(-0.0625, 0.0625),
                                          scale_limit=(-0.1, 0.1),
                                          rotate_limit=(-0, 0))
            #             x = randomHueSaturationValue(x,
            #                                hue_shift_limit=(-50, 50),
            #                                sat_shift_limit=(-5, 5),
            #                                val_shift_limit=(-15, 15))
            X_batch.append(x)
            Y_batch.append(y)
        X = np.asarray(X_batch, dtype=np.float32)
        Y = np.asarray(Y_batch, dtype=np.float32)
        #return X,Y
        yield X, Y

def delete_files_from_directory(dirpath):
    for filename in os.listdir(dirpath):
        filepath = os.path.join(dirpath, filename)
        try:
            shutil.rmtree(filepath)
        except OSError:
            os.remove(filepath)

def generate_training_batch_augmentor(data, batch_size):
    while True:
        X_batch = []
        Y_batch = []
        batch_ids = np.random.choice(data,size=batch_size,replace=False)

        delete_files_from_directory(temp_destination_Input)
        delete_files_from_directory(temp_destination_GT)
        for idx, img_id in enumerate(batch_ids):
            mask_id=img_id[:-3]+'png'
            shutil.copy(Training_Input+img_id,temp_destination_Input)
            shutil.copy(Training_GroundTruth+mask_id,temp_destination_GT)


        p = Aug.Pipeline(temp_destination_Input)

        p.ground_truth(temp_destination_GT)
        condition_thresh=np.random.random()

        # p.rotate(probability=0.5, max_left_rotation=10, max_right_rotation=10)
        # p.flip_left_right(probability=0.5)  # random flip operation
        # p.flip_top_bottom(probability=0.5)  # random flip operation
        # p.random_distortion(probability=0.5, grid_width=10, grid_height=10, magnitude=8)

        if condition_thresh < thr_pp:
             p.rotate(probability=0.8, max_left_rotation=10, max_right_rotation=10)
        elif 2*thr_pp>condition_thresh >=thr_pp:
             p.flip_left_right(probability=0.8)  # random flip operation
             p.flip_top_bottom(probability=0.8)  # random flip operation
        elif 3*thr_pp>condition_thresh >=2*thr_pp:
             p.random_distortion(probability=0.8, grid_width=8, grid_height=8, magnitude=5)
        elif 4*thr_pp>condition_thresh >=3*thr_pp:
             p.rotate(probability=0.5, max_left_rotation=10, max_right_rotation=10)
             p.zoom(probability=0.8, min_factor=1.1, max_factor=1.5)
        elif 5*thr_pp>condition_thresh >=4*thr_pp:
             p.rotate(probability=0.5, max_left_rotation=10, max_right_rotation=10)
             p.skew(probability=0.8)
        elif 6*thr_pp>condition_thresh >=5*thr_pp:
             p.rotate(probability=0.5, max_left_rotation=10, max_right_rotation=10)
             p.flip_left_right(probability=0.5)  # random flip operation
             p.flip_top_bottom(probability=0.5)  # random flip operation
        elif 7*thr_pp>condition_thresh >=6*thr_pp:
            p.random_distortion(probability=0.5, grid_width=8, grid_height=8,
                                magnitude=5)  # randome elastic transformation
            p.random_contrast(probability=0.8, min_factor=0.9, max_factor=1.4)
        else:
             p.rotate(probability=0.5, max_left_rotation=10, max_right_rotation=10)
             p.random_distortion(probability=0.8, grid_width=10, grid_height=10, magnitude=8)

        p.sample(batch_size*2)
        temp_img_ids=get_img_id_xu()
        for idx,temp_id in enumerate(temp_img_ids):
            x=imread(os.path.join(temp_input_augmented,temp_id),as_grey=False)
            y=imread(os.path.join(temp_input_augmented,'_groundtruth_(1)_temp_Input_'+temp_id[20:]),as_grey=False)

            # always keep the gt mask one-hot point encoding
            ind_max = np.argmax(y, axis=2)
            y_mask = np.zeros((HEIGHT, WIDTH, CHANNELS), np.float32)
            xx = np.linspace(0, HEIGHT - 1, HEIGHT, dtype=np.int)
            yy = np.linspace(0, WIDTH - 1, WIDTH, dtype=np.int)
            xv, yv = np.meshgrid(xx, yy)
            y_mask[yv, xv, ind_max] = 1

            X_batch.append(x)
            Y_batch.append(y_mask)

        X = np.asarray(X_batch, dtype=np.float32)
        Y = np.asarray(Y_batch, dtype=np.float32)

        #return X,Y
        yield X, Y

def generate_validation_batch_xu(data, batch_size):
    while True:
        X_batch = []
        Y_batch = []
        batch_ids = np.random.choice(data,
                                     size=batch_size,
                                     replace=False)
        for idx, img_id in enumerate(batch_ids):
            x=imread(os.path.join(Training_Input,img_id),as_grey=False)
            mask_id = img_id[:-3] + 'png'
            y = imread(os.path.join(Training_GroundTruth,mask_id),as_grey=False)

            ind_max = np.argmax(y, axis=2)
            y_mask = np.zeros((HEIGHT, WIDTH, CHANNELS), np.float32)
            xx = np.linspace(0, HEIGHT - 1, HEIGHT, dtype=np.int)
            yy = np.linspace(0, WIDTH - 1, WIDTH, dtype=np.int)
            xv, yv = np.meshgrid(xx, yy)
            y_mask[yv, xv, ind_max] = 1

            X_batch.append(x)
            Y_batch.append(y_mask)
        X = np.asarray(X_batch, dtype=np.float32)
        Y = np.asarray(Y_batch, dtype=np.float32)
        yield X, Y
        #return X, Y

def generate_validation_batch(data, batch_size,train_imgs,train_masks):
    while True:
        X_batch = []
        Y_batch = []
        batch_ids = np.random.choice(data,
                                     size=batch_size,
                                     replace=False)
        for idx, img_id in enumerate(batch_ids):
            #x = get_image(img_id)
            x=train_imgs[img_id]
            mask_id = img_id[:-3] + 'png'
            #y = get_mask(mask_id)
            y=train_masks[mask_id]
            X_batch.append(x)
            Y_batch.append(y)
        X = np.asarray(X_batch, dtype=np.float32)
        Y = np.asarray(Y_batch, dtype=np.float32)
        yield X, Y

def train_segmentation():
    img_ids, mask_ids = get_img_id()

    # train_imgs = {}
    # for img_id in tqdm(img_ids):  # tqdm show progress for loop
    #     train_imgs[img_id] = load_image_disk(img_id)
    #
    # train_masks = {}
    # for img_id in tqdm(mask_ids):
    #     train_masks[img_id] = load_mask_disk(img_id)




    # Training new model
    ts = str(int(time.time()))
    model_name = 'U-net-NS'
    steps_per_epoch = int(len(img_ids*10) / BATCH_SIZE)
    run_name = 'model={}-batch_size={}-num_epoch={}-steps_per_epoch={}-ts={}'.format(model_name,
                                                                                     BATCH_SIZE,
                                                                                     num_epochs,
                                                                                     steps_per_epoch,
                                                                                     ts)
    tensorboard_loc = os.path.join(TENSORBOARD_PATH, run_name)
    checkpoint_loc = os.path.join(MODELS_PATH, 'model-{}-weights.h5'.format(ts))

    earlyStopping = EarlyStopping(monitor='val_loss', # quantity to be monitored
                                  patience=10,       # number of epochs with no improvement after which training will be stopped
                                  verbose=1,        # decides what to print
                                  min_delta=0.00001, # minimum change to qualify as an improvement
                                  mode='min', )     # min mode: training will stop when the quantity monitored has stopped decreasing

    modelCheckpoint = ModelCheckpoint(checkpoint_loc,   # path to save the model
                                      monitor='val_loss', # quantity to monitor
                                      save_best_only=True, # if ture, the latest best model will not be overwritten
                                      mode='min',          # the decision to overwrite based on
                                      verbose=1,
                                      save_weights_only=True)

    tensorboard = TensorBoard(log_dir=tensorboard_loc, histogram_freq=0, write_graph=True, write_images=True)

    callbacks_list = [modelCheckpoint, earlyStopping, tensorboard]

    model = get_unet()

    print(model.summary())
    get_model_memory_usage(BATCH_SIZE, model)

    train_ids, validation_ids = model_selection.train_test_split(img_ids, random_state=48, test_size=0.2)
    train_generator = generate_training_batch_augmentor(train_ids, BATCH_SIZE)
    valid_generator = generate_validation_batch_xu(validation_ids, BATCH_SIZE)
    #train_generator = generate_training_batch(train_ids, BATCH_SIZE,train_imgs,train_masks)
    #valid_generator = generate_validation_batch(validation_ids, BATCH_SIZE,train_imgs,train_masks)
    VALIDATION_STEPS = int(len(validation_ids) / BATCH_SIZE)

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

if __name__ == '__main__':
    train_segmentation()