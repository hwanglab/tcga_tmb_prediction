import os
#import h5py
from keras.models import model_from_json
from keras.preprocessing import image
from skimage.io import imshow
import numpy as np
from scipy.io import loadmat
import glob

Test_data=['E:/Hongming/projects/tcga-bladder-mutationburden/deep_learning_data/data01/validation/low/',
           'E:/Hongming/projects/tcga-bladder-mutationburden/deep_learning_data/data01/validation/mid/',
           'E:/Hongming/projects/tcga-bladder-mutationburden/deep_learning_data/data01/validation/high/']

model_location='model_first_try.json'
model_weights='frist_try.h5'
patientDI=loadmat('E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/'+'blca_MutBurdens.mat')
cvdiv=loadmat('E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/'+'cvIndices.mat')
lbpPath='E:/Hongming/projects/tcga-bladder-mutationburden/featureoutput/10)AP/1)lbp/'
fold=[1]

def get_img_id():
    test_ids=[]
    for k in range(1,2):
        path_temp=Test_data[k]
        images=os.listdir(path_temp)
        for image_name in images:
            if '.png' or '.jpg' in image_name:
                test_ids.append(image_name)

    return test_ids

def load_model_prediction():
    json_file=open(model_location, 'r')
    model_json= json_file.read()
    json_file.close()
    model = model_from_json(model_json)
    model.load_weights(model_weights)

    cvdivff=cvdiv['cvIndices']
    pID=patientDI['blca_MutBurdens']
    tnum=0
    cnum=0
    for k in range(0,len(cvdivff)):
        temp=cvdivff[k]
        if temp==fold:
            pID_temp=pID[k,0] # patient ID
            cc_temp=pID[k,2]  # class type
            matff = glob.glob(lbpPath + pID_temp[0] + '*.mat')
            if len(matff) == 1:
                tnum += 1
                if cc_temp == 'Low':
                    index = tmb_prediction(model, pID_temp[0], Test_data[0],matff[0])
                    if index == 1:
                        cnum += 1
                elif cc_temp == 'Mid':
                    index = tmb_prediction(model, pID_temp[0], Test_data[1],matff[0])
                    if index == 2:
                        cnum += 1
                elif cc_temp == 'High':
                    index = tmb_prediction(model, pID_temp[0], Test_data[2],matff[0])
                    if index == 0:
                        cnum += 1
                else:
                    print('impossible!!!!')

    print('accuray is {}'.format(cnum/tnum))

def tmb_prediction(model,pID_temp,image_path,matfile):
    images = glob.glob(image_path + pID_temp + '*.png')
    lbpf=loadmat(matfile)
    feat_out=lbpf['feat_out']
    coeff=feat_out[:,-1]/sum(feat_out[:,-1])
    pred=[]
    for image_name in images:
        if pID_temp in image_name:
            img=image.load_img(image_name,target_size=(252,252))
            x=image.img_to_array(img)
            x=np.expand_dims(x,axis=0)
            x=x/255
            pred_lab=model.predict(x)
            pred.append(np.ndarray.tolist(pred_lab[0]))

    pred=np.asarray(pred)
    ss=pred.shape
    ww=np.tile(np.array([coeff]).transpose(), (1, ss[1]))
    predw=np.multiply(pred,ww)
    predf=np.sum(predw,axis=0)
    index=np.argmax(predf)

    return index

if __name__=='__main__':
    #test_ids=get_img_id()
    load_model_prediction()