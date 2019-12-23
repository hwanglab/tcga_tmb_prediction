import os
import numpy as np
from keras.models import load_model
import glob
import scipy.io as sio


test_path=['E:/Hongming/projects/tcga-bladder-mutationburden/featureoutput/2)mutation_prediction/testing/low/',
           'E:/Hongming/projects/tcga-bladder-mutationburden/featureoutput/2)mutation_prediction/testing/mid/',
           'E:/Hongming/projects/tcga-bladder-mutationburden/featureoutput/2)mutation_prediction/testing/high/']

MODELPATH ='E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/step4)_wsi_classification/Assets/models/'

if __name__=='__main__':
    print('-' * 30)
    print('Loading saved model...')
    print('-' * 30)

    #model = load_model(os.path.join(MODELPATH, 'model-1535553973.h5'),
    #                   custom_objects={'dice_coef_loss': dice_coef_loss})
    model = load_model(os.path.join(MODELPATH, 'model-1537372786.h5'))

    x_train = []
    y_train = []
    x_valid = []
    y_valid = []
    for ind in range(0, len(test_path)):
        path_ind = test_path[ind]
        matfile = glob.glob(path_ind + "*.mat")

        for k in range(0, len(matfile)):
            mat_contents = sio.loadmat(matfile[k])
            x_train.append(mat_contents['feat_out'])
            y_train.append(ind)

    X_Train = np.asarray(x_train)
    y_out = model.predict(X_Train)
    y=np.argmax(y_out,axis=1)

    t=0