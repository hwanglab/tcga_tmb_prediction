# tcga_tmb_prediction
Using transfer learning to predict tumor mutational burden from whole slide images for bladder cancer patients

Main author: Hongming Xu, CCF, 2019

Emails: xuh3@ccf.org or mxu@ualberta.ca

## Codes:
-- Matlab: the majority part of image analysis and classification

-- Python: use pre-trained models to extract tumor tile features

-- R: to plot KM survival curves

## Reqirements:
-- pathology slides: you need to get access tcga_blca .svs pathology slides from tcga data portal (freely downloading) OR you can get access them from hwang_lab shared disk (for hwang lab memembers)

-- matlab toolbox: you need to download matlab openslide toolbox (freely online)

-- python toolbox: you need to install tensorflow+keras (freely online)

## Usages:
-- step00: pre-processing during our project develompent, you do not need to run this folder for repeating our studies

-- step01: we trained svm classifier to detect tumor regions
