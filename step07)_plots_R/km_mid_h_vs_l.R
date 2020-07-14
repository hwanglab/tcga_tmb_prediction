# load packages

library(survival)
library(survminer)
library(dplyr)
library(readxl)
library(R.matlab)
library(openxlsx)

data_path0<-"E:/Hongming/projects/tcga-bladder-mutationburden/tcga_tmb_prediction/step06)_plot_km_curves/"
pid_pred<-read_excel(paste(data_path0,"pid_tmb_pred.xlsx",sep=""),sheet = "Sheet1")
pid<-pid_pred$`patient_names`
pred<-pid_pred$preds

blca_pred<-data.frame("patientID"=Reduce(rbind,pid),"label_class"=Reduce(rbind,pred))

data_path<-"E:/Hongming/projects/tcga-bladder-mutationburden/tcga_tmb_prediction/"
blca_data<-read_excel(paste(data_path,"Table_S1.2017_08_05.xlsx",sep=""),sheet="Master table")

futime<-vector()
fustat<-vector()
for (nn in 1:length(pid))
{
  temp_pID=as.character(pid[nn])
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


row_ind2<-which(blca_pred$futime<0)
if (length(row_ind2)>0){
  blca_pred<-blca_pred[-c(row_ind2),]
}

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
ggsurvplot(fit1,pval = TRUE,risk.table = TRUE, xlab="Time in months")