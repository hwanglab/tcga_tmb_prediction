import os
import numpy as np
from keras.models import load_model
import glob
import scipy.io as sio
from skimage.io import imsave, imread,imshow

test_path=['E:/Hongming/projects/tcga-bladder-mutationburden/tumor_patches_5x/testing/low/',
           'E:/Hongming/projects/tcga-bladder-mutationburden/tumor_patches_5x/testing/mid/',
           'E:/Hongming/projects/tcga-bladder-mutationburden/tumor_patches_5x/testing/high/']

MODELPATH ='E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/step4)_wsi_classification/Assets/vggmodels/'

if __name__=='__main__':
    print('-' * 30)
    print('Loading saved model...')
    print('-' * 30)

    #model = load_model(os.path.join(MODELPATH, 'model-1535553973.h5'),
    #                   custom_objects={'dice_coef_loss': dice_coef_loss})
    model = load_model(os.path.join(MODELPATH, 'model-1537481410.h5'))


    y_train = []
    x_valid = []
    y_valid = []
    for ind in range(0, len(test_path)):
        path_ind = test_path[ind]
        jpgfile = glob.glob(path_ind + "*.jpg")

        x_test = []
        for k in range(0, len(jpgfile)):
            x = imread(os.path.join(path_ind, jpgfile[k]), as_grey=False)
            x_test.append(x)


        X_Test = np.asarray(x_test)
        y_out = model.predict(X_Test)
    #y=np.argmax(y_out,axis=1)

        t=0