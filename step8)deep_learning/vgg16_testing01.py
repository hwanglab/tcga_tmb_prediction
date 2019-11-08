import os
#import h5py
from keras.models import model_from_json
from keras.preprocessing import image
from skimage.io import imshow
import numpy as np

#Test_data=['E:/Hongming/deep_learning/keras_learning/data/test/cats/',
#           'E:/Hongming/deep_learning/keras_learning/data/test/dogs/',
#           'E:/Hongming/deep_learning/keras_learning/data/test/flowers/']

Test_data=['E:/Hongming/projects/tcga-bladder-mutationburden/deep_learning_data/data01/validation/low/',
           'E:/Hongming/projects/tcga-bladder-mutationburden/deep_learning_data/data01/validation/mid/',
           'E:/Hongming/projects/tcga-bladder-mutationburden/deep_learning_data/data01/validation/high/']

model_location='model_first_try.json'
model_weights='frist_try.h5'

def get_img_id():
    test_ids=[]
    for k in range(1,2):
        path_temp=Test_data[k]
        images=os.listdir(path_temp)
        for image_name in images:
            if '.png' or '.jpg' in image_name:
                test_ids.append(image_name)

    return test_ids

def load_model_prediction(test_ids):
    json_file=open(model_location, 'r')
    model_json= json_file.read()
    json_file.close()
    model = model_from_json(model_json)
    model.load_weights(model_weights)

    for k in range(0,len(test_ids)):
        temp=Test_data[1]+test_ids[k]
        img=image.load_img(temp,target_size=(252,252))
        x=image.img_to_array(img)  # converts a PIL image instance to a numpy array, returns a 3D numpy array
        x=np.expand_dims(x,axis=0)
        x=x/255

        pred_lab=model.predict(x)
        print(pred_lab)

if __name__=='__main__':
    test_ids=get_img_id()
    load_model_prediction(test_ids)