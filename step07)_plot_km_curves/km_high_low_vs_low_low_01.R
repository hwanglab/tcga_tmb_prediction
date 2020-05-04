# function to plot KM curves
# high_low vs low_low
# without considering intermediate level patients



# load packages
library(survival)
library(survminer)
library(dplyr)
library(readxl)
library(R.matlab)
library(openxlsx)

data_path0<-"E:/Hongming/projects/tcga-bladder-mutationburden/heatmap_blca/"
pid_pred_gt<-read_excel(paste(data_path0,"pid_pred_gt.xlsx",sep=""),sheet = "Sheet1")


path <- system.file("mat-files", package = "R.matlab")
pathname <- file.path("E:/Hongming/projects/tcga-bladder-mutationburden/tcga_tmb_prediction/step06)_plot_km_curves/features", "feat.mat")
blca_pre<-readMat(pathname)
pid<-vector()
plabel<-vector()
pvalue<-vector()
plabel0<-vector()
plabel1<-vector()
plabel<-vector()

for (kk in 1:length(blca_pre[[1]])/3){
  pid[kk]<-blca_pre[[1]][[kk]]
  pvalue[kk]<-blca_pre[[1]][[kk+length(blca_pre[[1]])/3]]
  plabel0[kk]<-blca_pre[[1]][[kk+length(blca_pre[[1]])*2/3]]
    
}

#tt<-median(as.numeric(pvalue))
tt<-quantile(as.numeric(pvalue),c(0.3,0.4,0.5,0.6,0.7))
plabel1<-(pvalue>tt[3])
plabel1[plabel1==TRUE]<-'High'
plabel1[plabel1==FALSE]<-'Low'

pwex0<-vector()
for (kk in 1:length(pid)){
  temp_pID<-substring(as.character(pid[kk]),1,23)
  for (jj in 1:length(pid_pred_gt$`Patient ID`))
    if (temp_pID==substring(as.character(pid_pred_gt$`Patient ID`[jj]),2,24)){
      plabel[kk]<-paste(pid_pred_gt$Automatic_Predictions[jj],'-',plabel1[kk],sep="")
      pwex0[kk]<-pid_pred_gt$Ground_Truths[jj]
      break
    }
}


# pid2<-vector()
# plabel2<-vector()
# pwex<-vector()
# ind=1
# for (kk in 1:length(pid)){
#   if (plabel[kk]=="'High'-Low") {
#     pid2[ind]<-substring(pid[kk],1,12)
#     plabel2[ind]<-plabel[kk]
#     pwex[ind]<-pwex0[kk]
#     ind<-ind+1
#   }
# }
# blca<-data.frame("patientID"=Reduce(rbind,pid2),"label_class"=Reduce(rbind,plabel2),"wex_class"=Reduce(rbind,pwex))
# write.xlsx(blca,"blca.xlsx")



pid2<-vector()
plabel2<-vector()
ind=1
for (kk in 1:length(pid)){
  if (plabel[kk]=="'High'-Low" | plabel[kk]=="'Low'-Low") {
    pid2[ind]<-pid[kk]
    plabel2[ind]<-plabel[kk]
    ind<-ind+1
  }
}

# pid2<-vector()
# plabel2<-vector()
# ind=1
# for (kk in 1:length(pid)){
#   if (plabel[kk]=="'High'-High" | plabel[kk]=="'Low'-Low") {
#     pid2[ind]<-pid[kk]
#     plabel2[ind]<-plabel[kk]
#     ind<-ind+1
#   }
# }

blca_pred<-data.frame("patientID"=Reduce(rbind,pid2),"label_class"=Reduce(rbind,plabel2))
#write.xlsx(blca,"blca.xlsx")


data_path<-"E:/Hongming/projects/tcga-bladder-mutationburden/tcga_tmb_prediction/"
blca_data<-read_excel(paste(data_path,"Table_S1.2017_08_05.xlsx",sep=""),sheet="Master table")


futime<-vector()
fustat<-vector()
for (nn in 1:length(blca_pred$patientID))
{
  temp_pID=as.character(blca_pred$patientID[nn])
  for (kk in 1:length(blca_data$`Case ID`))
  {
    if (substring(temp_pID,1,12)==as.character(blca_data$`Case ID`[kk]))
    {
      futime<-c(futime,as.numeric(blca_data$`Combined days to last followup or death`[kk])/30)
      fustat<-c(fustat,as.character(blca_data$`Vital status`[kk]))
      break
    }
  }
}

blca_pred$futime<-futime
blca_pred$fustat<-fustat

#row_ind<-which(is.na(blca_pred$futime))
#blca_pred<-blca_pred[-c(row_ind),]   # remove two not available patients

row_ind2<-which(blca_pred$futime<0)
if (length(row_ind2)>0){
  blca_pred<-blca_pred[-c(row_ind2),]
}

#rownames(blca2)<-NULL # reorder row number

# computer median survival time for different group of patients
class1_surv<-c()
class2_surv<-c()
for (m in 1:length(blca_pred$patientID))
{
  temp=as.character(blca_pred$label_class[m])
  if (temp=="'High'-Low"){
    class1_surv<-c(class1_surv,blca_pred$futime[m])
  } else if(temp=="'Low'-Low"){
    class2_surv<-c(class2_surv,blca_pred$futime[m])
  }else {print("not possible!!!!")}
}
surv1_mean<-mean(class1_surv)
surv1_median<-median(class1_surv)

surv2_mean<-mean(class2_surv)
surv2_median<-median(class2_surv)

# dichotomize the dead/alive and change data labels
# dead: censored 1, alive 0
blca_pred$fustat <- factor(blca_pred$fustat,
                       levels=c("Dead","Alive"),
                       labels=c("1","0"))

blca_pred$futime<-as.numeric(as.character(blca_pred$futime)) # note that: must be numeric type
blca_pred$fustat<-as.numeric(as.character(blca_pred$fustat)) # note that: must be numeric type


# fit survival data using the kaplan-Meier method
surv_object<-Surv(time=blca_pred$futime,event=blca_pred$fustat)
# surv_object

fit1<-survfit(surv_object~label_class,data=blca_pred)
# summary(fit1)
# 

#setEPS()
#postscript("whatever.eps")

ggsurvplot(fit1,pval = TRUE,
           #risk.table = TRUE,
           legend=c(0.8,0.9),
           legend.labs=c("High-Low (61)","Low-Low (63)"),
           legend.title="Categories",
           xlab='Time in Months')+ggtitle("TCGA Bladder Cohort")

# univariate cox analysis for a certain category
res.cox<-coxph(surv_object~label_class,data=blca_pred)
res<-summary(res.cox)
#dev.off()





