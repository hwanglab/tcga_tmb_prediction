# tcga_tmb_prediction
Deep learning approach to predict tumor mutation burden (TMB) and delineate its spatial heterogeneity from whole slide images.

Main author: Hongming Xu, CCF, 2019

Emails: xuh3@ccf.org or mxu@ualberta.ca

## Codes:
-- Matlab: the majority part of image analysis and classification

-- Python: use pre-trained models to extract tumor tile features

-- R: plot KM survival curves

## Reqirements:
-- pathology slides: you need to get access tcga_blca .svs pathology slides from tcga data portal (freely downloading) OR you can get access them from hwang_lab shared disk (for hwang lab members)

-- matlab toolbox: you need to download matlab openslide toolbox (freely online)

-- python toolbox: you need to install tensorflow+keras (freely online)

## Simple Example:

-- You can get the basic ideas for our study by running codes in: example_tmb_prediction folder

## Explanations of different folders:

-- step00: contains the code to divide patients into low, mid and high TMB groups. It is used as the preprocessing, you do not need to run this folder for repeating our studies

-- step01: for simplicity, we trained a svm classifier to detect tumor regions (which is used in the example). You can skip this folder to run other steps too. (deep learning based tumor detector codes are also included inside this folder)

-- step02: input the .svs whole slides images, we detect and save selected tumor tiles by ap clustering

-- step03: extract feaures from selected tumor tiles

-- step04: tmb prediction classifications

-- step05: generate heatmaps and compute entropy

-- step06: generate tmb prediction roc curves

-- step07: generate km curves

## Notes:
-- utility_funcs：it saves all reqiured matlab functions for this study


