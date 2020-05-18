# function to plot KM curves
# high vs low for all patients

# load packages
library(survival)
library(survminer)
library(dplyr)
library(readxl)
library(R.matlab)
library(openxlsx)

data_path0<-"E:/Hongming/projects/tcga-bladder-mutationburden/tcga_tmb_prediction/step07)_plot_km_curves/patients_all/"
pid_tmb_pred<-read_excel(paste(data_path0,"pid_tmb_pred.xlsx",sep=""),sheet = "Sheet1")

#data_path0<-"E:/Hongming/projects/tcga-bladder-mutationburden/heatmap_blca/"
#pid_pred_gt<-read_excel(paste(data_path0,"pid_pred_gt.xlsx",sep=""),sheet = "Sheet1")

data_path1<-"E:/Hongming/projects/tcga-bladder-mutationburden/tcga_tmb_prediction/step07)_plot_km_curves/patients_all/"
entropy_all<-read_excel(paste(data_path1,"entropy_all.xlsx",sep=""),sheet = "Sheet1")


pvalue<-entropy_all$entropy
pid<-entropy_all$`patient ID`


plabel<-vector()
plabel1<-vector()



for (kk in 1:length(pid)){
  temp_pID<-substring(as.character(pid[kk]),2,24)
  for (jj in 1:length(pid_tmb_pred$`patient_names`))
    if (temp_pID==substring(as.character(pid_tmb_pred$`patient_names`[jj]),1,23)){
      plabel[kk]<-paste(pid_tmb_pred$preds[jj])
      break
    }
}

pid2<-vector()
plabel2<-vector()
ind=1
for (kk in 1:length(pid)){
  if (plabel[kk]=="High") {
    pid2[ind]<-substring(pid[kk],2,24)
    plabel2[ind]<-plabel[kk]
    ind<-ind+1
  } else if (plabel[kk]=="Low"){
    pid2[ind]<-substring(pid[kk],2,24)
    plabel2[ind]<-plabel[kk]
    ind<-ind+1
  } 
  else {
    print('skip!')
  }
}

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
      futime<-c(futime,as.numeric(blca_data$`Combined days to last followup or death`[kk])/30.0)
      fustat<-c(fustat,as.character(blca_data$`Vital status`[kk]))
      break
    }
  }
}

blca_pred$futime<-futime
blca_pred$fustat<-fustat

row_ind<-which(is.na(blca_pred$futime))
if (length(row_ind)>0) {
  blca_pred<-blca_pred[-c(row_ind),]   # remove two not available patients
}

row_ind2<-which(blca_pred$futime<0)
if (length(row_ind2)>0){
  blca_pred<-blca_pred[-c(row_ind2),]
}


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

#setEPS()
#postscript("whatever.eps")

ggsurvplot(fit1,pval = TRUE,
           risk.table = TRUE, 
           legend=c(0.8,0.9),
           #legend.labs=c("High-Low (76)","Low-Low (108)"),
           legend.title="Categories",
           xlab="Time in months")+ggtitle("Whole Bladder Cohort")

#dev.off()
ggsave(filename = "km_2class_all2.eps",
       fallback_resolution = 600,
       device = cairo_ps)



