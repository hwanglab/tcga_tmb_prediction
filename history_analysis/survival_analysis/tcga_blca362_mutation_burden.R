# load packages

library(survival)
library(survminer)
library(dplyr)
library(readxl)
library(R.matlab)
library(openxlsx)

data_path<-"../../"
blca_data<-read_excel(paste(data_path,"Table_S1.2017_08_05.xlsx",sep=""),sheet="Master table")
tmb<-blca_data$`total number Single Nucleotide Variants`
LI<-0.33
IH<-0.66
tLI<-quantile(tmb,probs=LI,names=FALSE)
tIH<-quantile(tmb,probs=IH,names=FALSE)


TMB_class<-cut(blca_data$`total number Single Nucleotide Variants`,
               breaks=c(-1,tLI,tIH,Inf),labels=c("Low","Mid","High"))

blca<-data.frame("patientID"=blca_data$`Case ID`,
                 "futime"=blca_data$`Combined days to last followup or death`,
                 "fustat"=blca_data$`Vital status`,
                 "TMB"=blca_data$`total number Single Nucleotide Variants`,
                 "TMB_class"=TMB_class)

#write.xlsx(blca,"blca.xlsx")

# dichotomize the dead/alive and change data labels
# dead: censored 1, alive 0
blca$fustat <- factor(blca$fustat,
                      levels=c("Dead","Alive"),
                      labels=c("1","0"))

image_path<-"E:/data/blca_mutationBurden/blca_wsi/"
files<-list.files(path=image_path,pattern = "\\.svs$")

selc<-rep(0,length(blca$patientID))

for (nn in 1:length(files)){
  pID_ss=substring(files[nn],1,12)
  for (mm in 1:length(blca$patientID)){
    pID_oo=as.vector(blca$patientID[mm])
    if (identical(pID_ss,pID_oo)){
      selc[mm]=1
      break
    }
  }
}

row_nos<-which(selc==0)
blca2<-blca[-c(row_nos),]   # remove not selected patients

futime_d<-as.numeric(as.character(blca2$futime))
row_ind<-which(is.na(futime_d))
blca2<-blca2[-c(row_ind),]   # remove two not available patients
blca2$futime<-as.numeric(as.character(blca2$futime))
blca2$fustat<-as.numeric(as.character(blca2$fustat))

# fit survival data using the kaplan-Meier method
surv_object<-Surv(time=blca2$futime,event=blca2$fustat)
# surv_object
# 
fit1<-survfit(surv_object~TMB_class,data=blca2)
# summary(fit1)
# 
ggsurvplot(fit1,data=blca2,pval=TRUE)

