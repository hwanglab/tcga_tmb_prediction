# load packages

library(survival)
library(survminer)
library(dplyr)
library(readxl)
library(R.matlab)
library(openxlsx)



path <- system.file("mat-files", package = "R.matlab")
pathname <- file.path("E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/step7)_plot_figures/", "xception_pca100.mat")
blca_pre<-readMat(pathname)
pid<-vector()
plabel<-vector()
for (kk in 1:length(blca_pre[[1]])/3){
  pid[kk]<-blca_pre[[1]][[kk]]
  plabel[kk]<-blca_pre[[1]][[kk+length(blca_pre[[1]])/3]]
  #plabel[kk]<-blca_pre[[1]][[kk+length(blca_pre[[1]])*2/3]]
    
}
blca_pred<-data.frame("patientID"=Reduce(rbind,pid),"TMB_class"=Reduce(rbind,plabel))
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
      futime<-c(futime,as.numeric(blca_data$`Combined days to last followup or death`[kk]))
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
blca_pred<-blca_pred[-c(row_ind2),]
#rownames(blca2)<-NULL # reorder row number

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

fit1<-survfit(surv_object~TMB_class,data=blca_pred)
# summary(fit1)
# 
ggsurvplot(fit1,
           legend=c(0.8,0.9),
           legend.title="Predictions by SVGG16",
           legend.labs=c("High TMB","Low TMB"),
           xscale="d_m",
           #risk.table = TRUE,
           #surv.median.line = "hv", # Add medians survival
           xlab="Time in months",
           pval=TRUE,
           #conf.int = TRUE,
           tables.theme = theme_cleantable(),
           palette = c("#FF9E29", "#86AA00"),
           #ggtheme = theme_gray()
           ) # Change ggplot2 theme


#ggsave(filename = "survival_svgg16.eps",
#       fallback_resolution = 600,
#       device = cairo_ps)




# # import the ovarian cancer dataset
# data(ovarian)
# glimpse(ovarian)
# 
# # dichotomize the age and change data labels
# ovarian$rx <- factor(ovarian$rx,
#                      levels=c("1","2"),
#                      labels=c("A","B"))
# 
# ovarian$resid.ds<-factor(ovarian$resid.ds,
#                          level=c("1","2"),
#                          labels=c("no","yes"))
# 
# ovarian$ecog.ps<-factor(ovarian$ecog.ps,
#                         levels=c("1","2"),
#                         labels=c("good","bad"))
# 
# # Data seems to be biomodal
# hist(ovarian$age)
# 
# ovarian<-mutate(ovarian,age_group=ifelse(age>=50,"old","young"))
# #ovarian <- ovarian %>% mutate(age_group = ifelse(age >=50, "old", "young"))
# # %>% is called multiple times to "chain" functions together
# # iris %>% head() %>% summary() is the same as: summary(head(iris))
# ovarian$age_group<-factor(ovarian$age_group)
# 
# # fit survival data using the kaplan-Meier method
# surv_object<-Surv(time=ovarian$futime,event=ovarian$fustat)
# surv_object
# 
# fit1<-survfit(surv_object~rx,data=ovarian)
# summary(fit1)
# 
# ggsurvplot(fit1,data=ovarian,pval=TRUE)