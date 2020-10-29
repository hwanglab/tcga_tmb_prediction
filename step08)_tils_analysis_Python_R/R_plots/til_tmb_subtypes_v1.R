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

data_path1<-"E:/Hongming/papers/2020 cancer patient tmb prediction/TMB_Manuscript/supplement/"
sup_tab<-read_excel(paste(data_path1,"temp1.xlsx",sep=""),sheet="Sheet1")
#til_tmbh_to_tum<-pid_tils$til_tmbh_to_tum

tils_tum=TRUE
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
  

  
  # plot the histogram
  #hist(as.numeric(unlist(pvalue)), main='Histogram of TILs densities',xlab="TILs Density", ylab="Frequency",breaks = 50)
  #abline(v=tt[2],col='red',lwd=3, lty=2)
  
  #plabel1<-as.vector(cut(as.numeric(unlist(pvalue)),breaks=c(-1,tt[1],tt[3],Inf),labels=c("Low","Mid","High")))
  
  plabel<-as.vector(df$'Fig5..a.')
  for (kk in 1:length(plabel)){
    if (plabel[kk]=='1' && plabel0[kk]=='High'){
      plabel0[kk]<-'HHL'
    }else{
      plabel0[kk]<-'Others'
    }
  }
  
  # generate sup files
  tils<-vector()
  figc<-vector()
  for (nn in 1 :length(sup_tab$`Case ID`)){
    temp_pID=as.character(sup_tab$`Case ID`[nn])
    
    ind<-which(pid==temp_pID)
    
    if (length(ind)==1) {
      tils<-c(tils,as.numeric(unlist(pvalue)[ind]))
      figc<-c(figc,plabel0[ind])
    } else {
      tils<-c(tils,'NA')
      figc<-c(figc,'NA')
    }
  }
  sup_tab['Tils_density']<-tils
  sup_tab['Fig5.(c)']<-figc
  #write.xlsx(sup_tab, paste(data_path1,"temp1.xlsx",sep=""), sheetName = "Sheet1")
  
  
  blca_pred<-data.frame("patientID"=Reduce(rbind,pid),"label_class0"=Reduce(rbind,plabel0))
  
  
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
  
  ggsurvplot(fit1,pval = TRUE,
                    #risk.table = TRUE, 
                    legend=c(0.75,0.85),
                    legend.labs=c("High-High-Low (HHL) (44)","Others (324)"),
                    legend.title="Categories",
                    xlab="Time in months")+ggtitle("TCGA BLCA (n=368)")
  #ggsave(file = "./2class_til_tum_tmb.pdf", print(survp,newpage = FALSE))
  ggsave(filename = "./km03.eps",
         fallback_resolution = 600,
         device = cairo_ps)
  
  
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
  ggsave(file = "./2class_til_tmbh_tumor_tmb.pdf", print(survp,newpage = FALSE))
  
}




