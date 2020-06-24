'''
By using this function, the blca wsi is divided into tiles, and each tile is predicted as tumor or non-tumor

Author: Hongming Xu, xuh3@ccf.org
'''

import numpy as np
import openslide
import scipy
import time
from TilingSchedule import TilingSchedule
from convert_schedule import convert_schedule
from wsi_preprocess_tissue import wsi_preprocess_tissue
#from overlap_heatmap import overlap_heatmap
import matplotlib.pyplot as plt
from PIL import Image


def wsi_tiling_prediction(Model, File, Magnification, Tile_size, MappingMag=1.25, Coverage=0.1,
                             zscore=False,scaling=False):
    """generating tiles and tumor prediction

    Parameters
    ----------
    Model: trained model
    File : str
        path and filename of slide.
    Magnification : double
        Desired magnification for tumor prediction.
    Tile : int
        Tile size used in disired magnification
    MappingMag: double
        low resolution magnification. Default value = 1.25.
    Coverage: double
        minimum percent of tile covered by tissue to be included.
        Ranges between [0,1). Default value = 0.1.

    Returns
    -------
    pmask : prediction result mask

    """

    # open image
    Slide = openslide.OpenSlide(File)

    # generate tiling schedule for desired sampling magnification
    Schedule = TilingSchedule(File, Magnification, Tile_size)

    # convert tiling schedule to low-resolution for tissue mapping
    lrSchedule = convert_schedule(Schedule, MappingMag)

    # get width, height of image at low-res reading magnification
    lrHeight = Slide.level_dimensions[lrSchedule.Level][1]
    lrWidth = Slide.level_dimensions[lrSchedule.Level][0]

    # read in whole slide at low magnification
    LR = Slide.read_region((0, 0), lrSchedule.Level, (lrWidth, lrHeight))

    # convert to numpy array and strip alpha channel
    LR = np.asarray(LR)
    LR = LR[:, :, :3]

    # resize if desired magnification is not provided by the file
    if lrSchedule.Factor != 1.0:
        LR = scipy.misc.imresize(LR, lrSchedule.Factor)


    # generate tissue foreground mask
    LRMask=wsi_preprocess_tissue(LR,210)

    pmask=np.zeros(LRMask.shape,dtype=np.float32)

    #prediction each foreground tile
    for i in range(Schedule.X.shape[0]-1):
        for j in range(Schedule.X.shape[1]-1):
            lrTileMask = LRMask[
                         i * lrSchedule.Tout:(i + 1) * lrSchedule.Tout,
                         j * lrSchedule.Tout:(j + 1) * lrSchedule.Tout].astype(np.uint8)
            TissueCount = sum(lrTileMask.flatten().astype(np.float)) / (lrSchedule.Tout ** 2)
            if TissueCount > Coverage:

                Tile = Slide.read_region((int(Schedule.X[i, j]),int(Schedule.Y[i, j])),Schedule.Level,(Schedule.Tout, Schedule.Tout))

                Tile = np.asarray(Tile)
                Tile = Tile[:, :, :3]

                # resize if desired magnification is not provided by the file
                if Schedule.Factor != 1.0:
                    Tile = scipy.misc.imresize(Tile, Schedule.Factor, interp='nearest').astype(np.float32)
                else:
                    Tile = Tile.astype(np.float32)

                if scaling==True:
                    Tile = scipy.misc.imresize(Tile, (224,224), interp='nearest').astype(float)

                if zscore==True:
                    # preprocessing
                    Tile*=1/255.0
                    Tile-=np.mean(Tile,keepdims=True)
                    Tile/=(np.std(Tile,keepdims=True)+1e-6)

                Tile = np.expand_dims(Tile, axis=0)

                pred_lab=Model.predict(Tile)

                pmask[i * lrSchedule.Tout:(i + 1) * lrSchedule.Tout,
                   j * lrSchedule.Tout:(j + 1) * lrSchedule.Tout]=pred_lab[0,0]

    return pmask

