# function to plot KM curves
# low_tmb_high_entropy vs low_tmb_low_entropy (low tmb is the ground truth)
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
pathname <- file.path("E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/step11)_tmb_map_process/features/", "feat.mat")
blca_pre<-readMat(pathname)
pid<-vector()
plabel<-vector()
pvalue<-vector()
#plabel0<-vector()
plabel1<-vector()


for (kk in 1:length(blca_pre[[1]])/3){
  pid[kk]<-blca_pre[[1]][[kk]]
  pvalue[kk]<-blca_pre[[1]][[kk+length(blca_pre[[1]])/3]]
  #plabel0[kk]<-blca_pre[[1]][[kk+length(blca_pre[[1]])*2/3]]
  
}


#tt<-median(as.numeric(pvalue))
tt<-quantile(as.numeric(pvalue),c(0.3,0.4,0.5,0.6,0.7))
plabel1<-(pvalue>tt[3])
plabel1[plabel1==TRUE]<-'High'
plabel1[plabel1==FALSE]<-'Low'

for (kk in 1:length(pid)){
  temp_pID<-substring(as.character(pid[kk]),1,23)
  for (jj in 1:length(pid_pred_gt$`Patient ID`))
    if (temp_pID==substring(as.character(pid_pred_gt$`Patient ID`[jj]),2,24)){
      plabel[kk]<-paste(pid_pred_gt$Ground_Truths[jj],'-',pid_pred_gt$Automatic_Predictions[jj],'-',plabel1[kk],sep="")
      break
    }
}


# pid2<-vector()
# plabel2<-vector()
# ind=1
# for (kk in 1:length(pid)){
#   if (plabel[kk]=="'High'-Low" | plabel[kk]=="'Low'-High") {
#     pid2[ind]<-pid[kk]
#     plabel2[ind]<-plabel[kk]
#     ind<-ind+1
#   }
# }


pid2<-vector()
plabel2<-vector()
ind=1
for (kk in 1:length(pid)){
  if (plabel[kk]=="'Low'-'High'-Low") {
    pid2[ind]<-pid[kk]
    plabel2[ind]<-plabel[kk]
    ind<-ind+1
  } else if(plabel[kk]=="'Low'-'Low'-Low"){
    pid2[ind]<-pid[kk]
    plabel2[ind]<-plabel[kk]
    ind<-ind+1
  } else if(plabel[kk]=="'Low'-'High'-High"){
    pid2[ind]<-pid[kk]
    plabel2[ind]<-plabel[kk]
    ind<-ind+1
  } 
  else {
    print('skip!')
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

data_path<-"E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/"
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

#row_ind<-which(is.na(blca_pred$futime))
#blca_pred<-blca_pred[-c(row_ind),]   # remove two not available patients

row_ind2<-which(blca_pred$futime<0)
if (length(row_ind2)>0){
  blca_pred<-blca_pred[-c(row_ind2),]
}

#rownames(blca2)<-NULL # reorder row number

# computer median survival time for different group of patients
#class1_surv<-c()
#class2_surv<-c()
#for (m in 1:length(blca_pred$patientID))
#{
#  temp=as.character(blca_pred$label_class[m])
#  if (temp=="'High'-Low"){
#    class1_surv<-c(class1_surv,blca_pred$futime[m])
#  } else if(temp=="'Low'-Low"){
#    class2_surv<-c(class2_surv,blca_pred$futime[m])
#  }else {print("not possible!!!!")}
#}
#surv1_mean<-mean(class1_surv)
#surv1_median<-median(class1_surv)

#surv2_mean<-mean(class2_surv)
#surv2_median<-median(class2_surv)

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
           #risk.table = TRUE, 
           legend=c(0.8,0.9),
           legend.labs=c("Low-High-High (20)","Low-High-Low (19)","Low-Low-Low (45)"),
           legend.title="Prediction categories",
           xlab="Time in months")+ggtitle("TCGA Bladder Cohort")

#dev.off()

#legend=c(0.8,0.9),
#legend.title="WES TMB Levels",
#legend.labs=c("High (126)","Low (121)"),
#risk.table = TRUE, 
#palette = c("#FF9E29", "#86AA00"),
#xlab="Time in months")+ggtitle("TCGA Bladder Cohort")

# ggsurvplot(fit1,
#            #legend=c(0.8,0.9),
#            legend.title="stratify by TMB-entropy",
#            #legend.labs=c("High TMB","Low TMB"),
#            xscale="d_m",
#            #risk.table = TRUE,
#            #surv.median.line = "hv", # Add medians survival
#            xlab="Time in months",
#            pval=TRUE,
#            conf.int = TRUE,
#            tables.theme = theme_cleantable(),
#            #palette = c("#FF9E29", "#86AA00"),
#            ggtheme = theme_gray()
#            ) # Change ggplot2 theme


#ggsave(filename = "survival_svgg16.eps",
#       fallback_resolution = 600,
#       device = cairo_ps)


