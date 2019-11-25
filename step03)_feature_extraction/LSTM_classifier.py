
import os
from keras.models import Sequential
from keras.layers import LSTM, Dense, Dropout
from keras.utils import to_categorical
import numpy as np
import glob
import time
import scipy.io as sio
from keras.optimizers import SGD,Adam
from keras.callbacks import EarlyStopping, ModelCheckpoint, TensorBoard
from keras import regularizers

feat_path=['E:/Hongming/projects/tcga-bladder-mutationburden/featureoutput/2)mutation_prediction/low/',
           'E:/Hongming/projects/tcga-bladder-mutationburden/featureoutput/2)mutation_prediction/mid/',
           'E:/Hongming/projects/tcga-bladder-mutationburden/featureoutput/2)mutation_prediction/high/']

#feat_path=['E:/Hongming/projects/tcga-bladder-mutationburden/featureoutput/3)mutation_prediction_cellular_features/low/',
#           'E:/Hongming/projects/tcga-bladder-mutationburden/featureoutput/3)mutation_prediction_cellular_features/mid/',
#           'E:/Hongming/projects/tcga-bladder-mutationburden/featureoutput/3)mutation_prediction_cellular_features/high/']

TENSORBOARD_PATH='E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/step4)_wsi_classification/Assets/tensorboard/'
MODELS_PATH ='E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/step4)_wsi_classification/Assets/models/'

train_p=0.8
data_dim = 113
#data_dim=180
timesteps = 10
num_classes = 3


x_train = []
y_train = []
x_valid = []
y_valid = []
for ind in range(0,len(feat_path)):
    path_ind=feat_path[ind]
    matfile=glob.glob(path_ind+"*.mat")

    for k in range (0,len(matfile)):
        mat_contents=sio.loadmat(matfile[k])

        if k<len(matfile)*train_p:

            x_train.append(mat_contents['feat_out'])
            #x_train.append(mat_contents['feat_cell'])
            y_train.append(ind)

        else:


            x_valid.append(mat_contents['feat_out'])
            #x_valid.append(mat_contents['feat_cell'])
            y_valid.append(ind)


X_Train=np.asarray(x_train)
ss=X_Train.shape
X_mean=np.mean(X_Train,axis=1)
X_mean2=np.repeat(X_mean[:,np.newaxis,:],ss[1],axis=1)
X_std=np.std(X_Train,axis=1)
X_std2=np.repeat(X_std[:,np.newaxis,:],ss[1],axis=1)
X_Train2=(X_Train-X_mean2)/X_std2

Y_Train=np.asarray(y_train)
Y_Train2=to_categorical(Y_Train,num_classes)

X_Valid=np.asarray(x_valid)
ssv=X_Valid.shape
X_mean_valid0=np.mean(X_Valid,axis=1)
X_mean_valid=np.repeat(X_mean_valid0[:,np.newaxis,:],ssv[1],axis=1)
X_std_valid0=np.std(X_Valid,axis=1)
X_std_valid=np.repeat(X_std_valid0[:,np.newaxis,:],ssv[1],axis=1)
X_Valid2=(X_Valid-X_mean_valid)/X_std_valid

Y_Valid=np.asarray(y_valid)
Y_Valid2=to_categorical(Y_Valid,num_classes)





ts = str(int(time.time()))
model_name = 'LSTM'
run_name = 'model={}-ts={}'.format(model_name,ts)
tensorboard_loc = os.path.join(TENSORBOARD_PATH, run_name)
checkpoint_loc = os.path.join(MODELS_PATH, 'model-{}-weights.h5'.format(ts))
earlyStopping = EarlyStopping(monitor='val_loss',  # quantity to be monitored
                                  patience=20,
                                  # number of epochs with no improvement after which training will be stopped
                                  verbose=1,  # decides what to print
                                  min_delta=0.0,  # minimum change to qualify as an improvement
                                  mode='min', )  # min mode: training will stop when the quantity monitored has stopped decreasing

modelCheckpoint = ModelCheckpoint(checkpoint_loc,  # path to save the model
                                      monitor='val_loss',  # quantity to monitor
                                      save_best_only=True,  # if ture, the latest best model will not be overwritten
                                      mode='min',  # the decision to overwrite based on
                                      verbose=1,
                                      save_weights_only=True)

tensorboard = TensorBoard(log_dir=tensorboard_loc, histogram_freq=0, write_graph=True, write_images=True)

callbacks_list = [modelCheckpoint, earlyStopping, tensorboard]

# expected input data shape: (batch_size, timesteps, data_dim)
model = Sequential()
model.add(LSTM(32, return_sequences=True,activation='relu',
               input_shape=(timesteps, data_dim)))  # returns a sequence of vectors of dimension 32
model.add(Dropout(0.05))
model.add(LSTM(32, return_sequences=True,activation='relu'))  # returns a sequence of vectors of dimension 32
regularizers.l1_l2(l1=0.005, l2=0.005)
#model.add(Dropout(0.2))
model.add(LSTM(16,activation='relu'))  # return a single vector of dimension 32
regularizers.l1_l2(l1=0.005, l2=0.005)
model.add(Dense(num_classes, activation='softmax'))

sgd = SGD(lr=0.001, decay=1e-6, momentum=0.9, nesterov=True)
model.compile(loss='categorical_crossentropy',
              optimizer=sgd,
              #optimizer=Adam(lr=1e-5),
              metrics=['accuracy'])

# Generate dummy training data
#x_train = np.random.random((1000, timesteps, data_dim))
#y_train = np.random.random((1000, num_classes))

# Generate dummy validation data
#x_val = np.random.random((100, timesteps, data_dim))
#y_val = np.random.random((100, num_classes))

model.fit(X_Train2, Y_Train2,
          batch_size=32, epochs=200,
          callbacks=callbacks_list,
          validation_data=(X_Valid2, Y_Valid2))

model_path = os.path.join(MODELS_PATH, 'model-{}.h5'.format(ts))
model.save(model_path)