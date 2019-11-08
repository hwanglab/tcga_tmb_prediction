library(openxlsx)
library(readxl)

data_path<-"E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/"
blca_data<-read_excel(paste(data_path,"Table_S1.2017_08_05.xlsx",sep=""),sheet="Master table")

pID=blca_data$'Case ID'
pCC=blca_data$'mRNA cluster'

data_path0<-"E:/Hongming/projects/tcga-bladder-mutationburden/heatmap_blca/"
pid_pred_gt<-read_excel(paste(data_path0,"pid_pred_gt.xlsx",sep=""),sheet = "Sheet1")
pID2<-pid_pred_gt$`Patient ID`
pID3<-substring(pID2,2,13)
ap_tmb<-pid_pred_gt$Automatic_Predictions
gt_tmb<-pid_pred_gt$Ground_Truths

data_path1<-"E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/step11)_tmb_map_process/features/"
pid_entropy<-read_excel(paste(data_path1,"feat.xlsx",sep = ""),sheet="Sheet1")

tt<-quantile(pid_entropy$entropy,c(0.3,0.4,0.5,0.6,0.7))
plabel1<-(pid_entropy$entropy>tt[3])
plabel1[plabel1==TRUE]<-'High'
plabel1[plabel1==FALSE]<-'Low'

H_H_LI<- H_H_LP<- H_H_L<- H_H_BS <- H_H_N <- H_H_ND<-0

H_L_LI<- H_L_LP<- H_L_L<- H_L_BS <- H_L_N <- H_L_ND<-0

L_H_LI<- L_H_LP<- L_H_L<- L_H_BS <- L_H_N <- L_H_ND<-0

L_L_LI<- L_L_LP<- L_L_L<- L_L_BS <- L_L_N <- L_L_ND<-0

for (i in 1:length(pid_entropy$`patient ID`))
{
  temp_id=pid_entropy$`patient ID`[i]
  pp<-substring(as.character(temp_id),2,13)
  
  ind0<-which(pp==pID3)
  tmb_cc<-ap_tmb[ind0]
  
  entropy_cc<-plabel1[i]
  
  ind1<-which(pp==pID)
  subtype_cc<-pCC[ind1]
  
  if (tmb_cc=="'High'")
  {
    if (entropy_cc=='High' & subtype_cc=='Luminal_infiltrated') {
      H_H_LI<-H_H_LI+1
    }
    else if(entropy_cc=='High' & subtype_cc=='Luminal_papillary'){
      H_H_LP<-H_H_LP+1
    }
    else if (entropy_cc=='High' & subtype_cc=='Luminal'){
      H_H_L<-H_H_L+1
    }
    else if (entropy_cc=='High' & subtype_cc=='Basal_squamous'){
      H_H_BS<-H_H_BS+1
    }
    else if (entropy_cc=='High' & subtype_cc=='Neuronal'){
      H_H_N<-H_H_N+1
    }
    else if (entropy_cc=='High' & subtype_cc=='ND'){
      H_H_ND<-H_H_ND+1
    }
    else if (entropy_cc=='Low' & subtype_cc=='Luminal_infiltrated') {
      H_L_LI<-H_L_LI+1
    }
    else if(entropy_cc=='Low' & subtype_cc=='Luminal_papillary'){
      H_L_LP<-H_L_LP+1
    }
    else if (entropy_cc=='Low' & subtype_cc=='Luminal'){
      H_L_L<-H_L_L+1
    }
    else if (entropy_cc=='Low' & subtype_cc=='Basal_squamous'){
      H_L_BS<-H_L_BS+1
    }
    else if (entropy_cc=='Low' & subtype_cc=='Neuronal'){
      H_L_N<-H_L_N+1
    }
    else if (entropy_cc=='Low' & subtype_cc=='ND'){
      H_L_ND<-H_L_ND+1
    }
    else
    {print('no possible---------')}
      
  }
  else if (tmb_cc=="'Low'")
  {
    if (entropy_cc=='High' & subtype_cc=='Luminal_infiltrated') {
      L_H_LI<-L_H_LI+1
    }
    else if(entropy_cc=='High' & subtype_cc=='Luminal_papillary'){
      L_H_LP<-L_H_LP+1
    }
    else if (entropy_cc=='High' & subtype_cc=='Luminal'){
      L_H_L<-L_H_L+1
    }
    else if (entropy_cc=='High' & subtype_cc=='Basal_squamous'){
      L_H_BS<-L_H_BS+1
    }
    else if (entropy_cc=='High' & subtype_cc=='Neuronal'){
      L_H_N<-L_H_N+1
    }
    else if (entropy_cc=='High' & subtype_cc=='ND'){
      L_H_ND<-L_H_ND+1
    }
    else if (entropy_cc=='Low' & subtype_cc=='Luminal_infiltrated') {
      L_L_LI<-L_L_LI+1
    }
    else if(entropy_cc=='Low' & subtype_cc=='Luminal_papillary'){
      L_L_LP<-L_L_LP+1
    }
    else if (entropy_cc=='Low' & subtype_cc=='Luminal'){
      L_L_L<-L_L_L+1
    }
    else if (entropy_cc=='Low' & subtype_cc=='Basal_squamous'){
      L_L_BS<-L_L_BS+1
    }
    else if (entropy_cc=='Low' & subtype_cc=='Neuronal'){
      L_L_N<-L_L_N+1
    }
    else if (entropy_cc=='Low' & subtype_cc=='ND'){
      L_L_ND<-L_L_ND+1
    }
    else
    {print('no possiblexxxxxxxx')}
  }
  else
  {print('no possible!')}
  
}
ctable <- matrix(c(H_H_LI, H_H_LP, H_H_L, H_H_BS, H_H_N,
                   H_L_LI, H_L_LP, H_L_L, H_L_BS, H_L_N,
                   L_H_LI, L_H_LP, L_H_L, L_H_BS, L_H_N,
                   L_L_LI, L_L_LP, L_L_L, L_L_BS, L_L_N), 
                      nrow = 4, ncol=5, byrow = TRUE,
                      dimnames =
                        list(c("H-H", "H-L", "L-H","L-L"),c("LI", "LP","L","BS","N")))

fisher.test(ctable,hybrid=TRUE)

chisq.test(ctable) 