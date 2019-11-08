This folder contains the program and intermediate results to build tumor_classifier

The running procedures:
1) main_feats_extraction_from_gt.m (save features for tumor patches)
2) main_feats_organization_tumor_nontumor.m (organize extracted image features)
3) two options:
   a) write the program two train SVM base tumor_classifer
   b) open matlab->apps->classification learner->train SVM classifer using your organized features, 
      then save the tumor detection clssifer using main_save_exportedmodel

The SVM_cubic_model and SVM_FineGaussian_model are two our trained SVM tumor detection classifier
You can use them to make tumor detection directly