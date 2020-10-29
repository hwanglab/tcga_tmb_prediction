'''
select HHHL subtype maps
'''

import pandas as pd
import shutil
import glob
df=pd.read_excel('./HHHL_subtypes.xlsx')

map_path='../step05)_heatmap_entropy/heatmap_blca/til_tmb_maps/'
des_path='./HHHL_subtype/'

pids=df['patientID'].tolist()
labels=df['label_class0'].tolist()

for kk, pp in enumerate(pids):
    if labels[kk]=='HHHL':
        imgs=glob.glob(map_path+pp+'*_probability.png')
        shutil.copy(imgs[0],des_path)
