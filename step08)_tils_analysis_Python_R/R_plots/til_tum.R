# purpose: plot til_tum_ratio km curves

library(survival)
library(survminer)
library(dplyr)
library(readxl)
#library(R.matlab)
library(openxlsx)

data_path0<-"E:/Hongming/projects/tcga-bladder-mutationburden/tcga_tmb_prediction/step08)_tils_analysis_Python_R/"
pid_tils<-read.csv(file = paste(data_path0,"tcga_blca_labels_v2.csv", sep=""))

data_path<-"E:/Hongming/projects/tcga-bladder-mutationburden/tcga_tmb_prediction/step00)_preprocessing/"
blca_data<-read_excel(paste(data_path,"Table_S1.2017_08_05.xlsx",sep=""),sheet="Master table")



#til_tmbh_to_tum<-pid_tils$til_tmbh_to_tum

tils_tum=FALSE
#tils_tmbh_tum=FALSE

if (tils_tum==TRUE){
  col_name<-"til_to_tum"
}else{
  col_name<-"til_tmbh_to_tum"
}

plots_km<-function(pid_tils,col_name,blca_data){
  
  df<-pid_tils[complete.cases(pid_tils[ , col_name]),]
  
  pvalue<-df[col_name]
  pid<-as.vector(df$'Case.ID')
  
  plabel0<-vector()
  
  tt<-quantile(as.numeric(unlist(pvalue)),c(0.333,0.5,0.667))
  plabel0<-(pvalue>tt[2])
  plabel0[plabel0==TRUE]<-'High'
  plabel0[plabel0==FALSE]<-'Low'
  
  plabel1<-as.vector(cut(as.numeric(unlist(pvalue)),breaks=c(-1,tt[1],tt[3],Inf),labels=c("Low","Mid","High")))
  
  
  blca_pred<-data.frame("patientID"=Reduce(rbind,pid),"label_class0"=Reduce(rbind,plabel0),"label_class1"=Reduce(rbind,plabel1))
  
  
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
  
  blca_pred$fustat <- factor(blca_pred$fustat,
                             levels=c("Dead","Alive"),
                             labels=c("1","0"))
  
  blca_pred$futime<-as.numeric(as.character(blca_pred$futime)) # note that: must be numeric type
  blca_pred$fustat<-as.numeric(as.character(blca_pred$fustat)) # note that: must be numeric type
  
  
  return(blca_pred)
  
}

blca_pred<-plots_km(pid_tils,col_name,blca_data)


if (tils_tum==TRUE){
  # fit survival data using the kaplan-Meier method
  surv_object<-Surv(time=blca_pred$futime,event=blca_pred$fustat)
  
  fit1<-survfit(surv_object~label_class0,data=blca_pred)
  
  survp<-ggsurvplot(fit1,pval = TRUE,
                    risk.table = TRUE, 
                    legend=c(0.8,0.9),
                    #legend.labs=c("High-High (97)","High-Low (76)","Low-High (87)","Low-Low (108)"),
                    legend.title="TILs density inside tumor",
                    xlab="Time in months")+ggtitle("Whole Bladder Cohort")
  ggsave(file = "./2class_til_tumor.pdf", print(survp,newpage = FALSE))
  
  fit1<-survfit(surv_object~label_class1,data=blca_pred)
  
  survp<-ggsurvplot(fit1,pval = TRUE,
                    risk.table = TRUE, 
                    legend=c(0.8,0.9),
                    #legend.labs=c("High-High (97)","High-Low (76)","Low-High (87)","Low-Low (108)"),
                    legend.title="TILs density inside tumor",
                    xlab="Time in months")+ggtitle("Whole Bladder Cohort")
  ggsave(file = "./3class_til_tumor.pdf", print(survp,newpage = FALSE))
}else{
  # fit survival data using the kaplan-Meier method
  surv_object<-Surv(time=blca_pred$futime,event=blca_pred$fustat)
  
  fit1<-survfit(surv_object~label_class0,data=blca_pred)
  
  survp<-ggsurvplot(fit1,pval = TRUE,
                    risk.table = TRUE, 
                    legend=c(0.8,0.9),
                    #legend.labs=c("High-High (97)","High-Low (76)","Low-High (87)","Low-Low (108)"),
                    legend.title="TILs density within tmb high inside tumor",
                    xlab="Time in months")+ggtitle("Whole Bladder Cohort")
  ggsave(file = "./2class_til_tmbh_tumor.pdf", print(survp,newpage = FALSE))
  
  fit1<-survfit(surv_object~label_class1,data=blca_pred)
  
  survp<-ggsurvplot(fit1,pval = TRUE,
                    risk.table = TRUE, 
                    legend=c(0.8,0.9),
                    #legend.labs=c("High-High (97)","High-Low (76)","Low-High (87)","Low-Low (108)"),
                    legend.title="TILs density with tmb high inside tumor",
                    xlab="Time in months")+ggtitle("Whole Bladder Cohort")
  ggsave(file = "./3class_til_tmbh_tumor.pdf", print(survp,newpage = FALSE))
}




