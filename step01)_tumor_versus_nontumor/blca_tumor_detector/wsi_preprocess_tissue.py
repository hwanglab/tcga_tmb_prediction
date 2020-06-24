'''
purpose: select tissue regions by simple thresholding

author: Hongming Xu, CCF
email: mxu@ualberta.ca
'''

from skimage import measure, morphology
import scipy
import numpy as np
import matplotlib.pyplot as plt

def wsi_preprocess_tissue(rgb,thrWhite=210):
    # step 1: select interested regions
    gray = rgb2gray(rgb)
    mask = np.zeros(np.shape(gray), dtype=np.uint8)
    mask[gray < thrWhite] = 1
    #mask=gray<thrWhite
    mask_label=measure.label(mask,background=0)

    properties=measure.regionprops(mask_label)
    thrNoise=round(max([prop.area for prop in properties])/3)
    bwTissue=morphology.binary_closing(mask,morphology.disk(1))
    bwTissue=morphology.remove_small_objects(bwTissue,thrNoise)

    return bwTissue


def rgb2gray(rgb):
    r, g, b = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]
    gray = 0.2989 * r + 0.5870 * g + 0.1140 * b

    return gray