#This folder contains the program and intermediate results to build tumor_classifier

------------------------------------------------------------------------------------
Input: given annotated tumor and on-tumor regions in the WSIs, this step builds a tumor detector using only LBP texture features with SVM classifier
Output: SVM Tumor Detector
------------------------------------------------------------------------------------
To reprodcuce my results, follow the following steps:

1) main_feats_extraction_from_gt.m (save features for tumor patches)
2) main_feats_organization_tumor_nontumor.m (organize extracted image features)
3) two options you can select:
   a) write the program to train SVM based tumor_classifer
   b) open matlab->apps->classification learner->train SVM classifer using your organized features,e.g., FT in my example 
      then save the tumor detection clssifer using main_save_exportedmodel
      (for simplicity, in my testing I use the option b)

The SVM_cubic_model and SVM_FineGaussian_model are two our trained SVM tumor detection classifier
You can use them to make tumor detection directly

svm_tumor_detection_evaluation.m -> testing trained tumor detectors
roc_tumor_detection_blca.m -> plot ROC curve for tumor detection
