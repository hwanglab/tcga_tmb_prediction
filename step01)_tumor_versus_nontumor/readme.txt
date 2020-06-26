#This folder contains the program and intermediate results to build tumor_classifier

# in our evaluations, we trained two tumor detectors:
(1) The first is LBP+SVM for tumor detection (see below explanations about our experiments)
(2) The designed CNN for tumor detection

------------------------------------------------------------------------------------
Input: given annotated tumor and on-tumor regions in the WSIs, this step builds a tumor detector using only LBP texture features with SVM classifier
Output: SVM Tumor Detector
------------------------------------------------------------------------------------

To reprodcuce my results, follow the following steps:

<<<<<<< HEAD
1) a_feats_extraction_gt.m (save features for tumor patches)
=======
1) a_feats_extraction_from_gt.m (save features for tumor patches)
>>>>>>> 9032d1c209d8f7ea8363d09a8871dd1e7d66c273
2) b_feats_organization.m (organize extracted image features)
3) two options you can select:
   a) write the program to train SVM based tumor_classifer
   b) open matlab->apps->classification learner->train SVM classifer using your organized features,e.g., FT in my example 
<<<<<<< HEAD
      then save the tumor detection clssifer using c_save_exportedmodel
=======
      then save the tumor detection clssifer using c_save_exportedmodel.m
>>>>>>> 9032d1c209d8f7ea8363d09a8871dd1e7d66c273
      (for simplicity, in my testing I use the option b)

The SVM_cubic_model and SVM_FineGaussian_model are two our trained SVM tumor detection classifier
You can use them to make tumor detection directly

svm_tumor_detection_evaluation.m -> testing trained tumor detectors
roc_tumor_detection_blca.m -> plot ROC curve for tumor detection
