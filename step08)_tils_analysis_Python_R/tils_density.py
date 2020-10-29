'''
purpose: quantify til densitys inside tumor regions

(1) for tcga blca
(2) for tcga luad
'''

import scipy.io
import glob
import matplotlib.pyplot as plt
from skimage import transform
import numpy as np
import pandas as pd

tcga_blca=False
tcga_luad=True
if tcga_blca==True:
    tmb_pred_path=['../step05)_heatmap_entropy/heatmap_blca/mat_files/high_low/',
                   '../step05)_heatmap_entropy/heatmap_blca/mat_files/mid/']
    til_pred_path='Y:/projects/data/tcga_blca_slide/til_maps/'
    master_table=pd.read_excel('./tcga_blca_labels.xlsx')
    pid=master_table['Case ID'].tolist()
elif tcga_luad==True:
    tmb_pred_path=['../../../tcga-lung-mutationburden/tcga_luad_tmb_predicion/4)plot_figures/heatmaps/']
    til_pred_path='Y:/projects/data/tcga_luad_slide/til_maps/'
    master_table=pd.read_csv('../../../tcga-lung-mutationburden/tcga_luad_tmb_predicion/5)survival_analysis/tcga_luad_table_pred.csv')
    pid = master_table['patient_names'].tolist()
else:
    raise RuntimeError('undefined option!!!!!!!!!!!!')


tils_thr=0.5

til_tum=[float("NaN")]*len(pid)
til_tmbh=[float("NaN")]*len(pid)

for temp_path in tmb_pred_path:
    mats=glob.glob(temp_path+'*.mat')
    for k, m in enumerate(mats):
        tmb=scipy.io.loadmat(m)
        img_name=mats[k].split('\\')[-1].split('.')[0] # \\ for win10 system
        til_img=glob.glob(til_pred_path+img_name+'_gray.png')
        if len(til_img)==1:
            til_map=plt.imread(til_img[0])

            tmb_map=transform.resize(np.nan_to_num(tmb['tmb_map']),til_map.shape[0:2],order=1) # 0:nearest neighbor

            #temp = np.nan_to_num(tmb['tmb_map'])
            #temp1 = temp[temp > 0]
            #tr=np.min(temp1)

            tum_mask=(tmb_map>0)*1 # tumor regions
            til_mask=(til_map[:,:,0]>tils_thr)*1 # tils regions

            til_tum_ratio=np.sum(np.logical_and(til_mask,tum_mask))/np.sum(tum_mask)
            if tcga_blca==True:
                pp=img_name[0:12]
            elif tcga_luad==True:
                pp=img_name
            else:
                raise RuntimeError('undefined option ....')
            ind=pid.index(pp)
            til_tum[ind]=til_tum_ratio

            tmbh_map=(tmb_map>0.5)*1
            til_tmbh_ratio=np.sum(np.logical_and(til_mask,tmbh_map))/np.sum(tum_mask)
            til_tmbh[ind]=til_tmbh_ratio

master_table['til_to_tum']=til_tum
master_table['til_tmbh_to_tum']=til_tmbh

if tcga_blca==True:
    master_table.to_csv('./tcga_blca_labels_v2.csv')
elif tcga_luad==True:
    master_table.to_csv('../../../tcga-lung-mutationburden/tcga_luad_tmb_predicion/5)survival_analysis/tcga_luad_table_pred.csv')
else:
    raise RuntimeError('undefined option ....')


