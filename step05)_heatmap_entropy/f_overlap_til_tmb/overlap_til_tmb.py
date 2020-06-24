'''
overlap til and tmb maps together for tcga_blca
'''
from scipy.io import loadmat
import glob
import numpy as np
import matplotlib.pyplot as plt
import skimage.transform
import skimage.filters
from skimage.io import imsave
import pandas as pd

tmb_path=['../heatmap_blca/mat_files/high_low/',
          '../heatmap_blca/mat_files/mid/']

tumor_path='../heatmap_blca/tumor_maps/'

til_maps=['Y:/projects/data/tcga_blca_slide/til_maps/wsi/',
          'Y:/projects/data/tcga_blca_slide/til_maps/wsi2/']

output_path='../heatmap_blca/til_tmb_maps/'

patientID=[]
til_tumor_density=[]
for path in tmb_path:
    tmbs=glob.glob(path+'*.mat')
    for i,name in enumerate(tmbs):
        tmp = loadmat(name)
        tmb_map=tmp['tmb_map']
        tmb_map[np.isnan(tmb_map)]=0

        til_map_name=name.split('\\')[-1].split('.')[0]+'_gray.png'  # win10 system
        try:
            til_map=plt.imread(til_maps[0]+til_map_name)
        except:
            try:
                til_map = plt.imread(til_maps[1] + til_map_name)
            except:
                print('no til map for img=%' % name)

        tumor_map = plt.imread(tumor_path+name.split('\\')[-1].split('.')[0]+'.png')
        tumor_map = skimage.transform.resize(tumor_map, til_map.shape[0:2],
                                         order=1)
        tmb_map = skimage.transform.resize(tmb_map, til_map.shape[0:2],
                                         order=1)  # 0 nearest-neighbor, 1: bi-linear (default)

        tumor_map=skimage.filters.gaussian(tumor_map,sigma=1.5)
        tmb_map=skimage.filters.gaussian(tmb_map,sigma=1.5)

        # save figure 1: probability maps
        til_map[:,:,1]=tmb_map
        til_map[:,:,2]=tumor_map
        #til_map[:,:,2]=np.zeros(tmb_map.shape)

        imsave(output_path + name.split('\\')[-1].split('.')[0] + '_probability.png', til_map)

        # save figure2: binary maps
        til_map2=til_map
        til_map2[:,:,0]=(til_map2[:,:,0]>=0.5)
        til_map2[:, :, 1] = (til_map2[:, :, 1] >= 0.5)
        til_map2[:, :, 2] = (til_map2[:, :, 2] >= 0.5)
        imsave(output_path + name.split('\\')[-1].split('.')[0] + '_binary.png', til_map2)

        tumor_mask=(til_map2[:,:,2]>=0.5)
        til_mask=(til_map2[:,:,0]>=0.5)

        til_tumor_density.append(np.sum(np.logical_and(til_mask,tumor_mask)) / np.sum(tumor_mask))
        patientID.append(name.split('\\')[-1].split('.')[0])

# save excel file
data={'img_id':patientID,'til_density':til_tumor_density}
df=pd.DataFrame(data)
df.to_csv(output_path+'til_density_of_tumor.csv')


