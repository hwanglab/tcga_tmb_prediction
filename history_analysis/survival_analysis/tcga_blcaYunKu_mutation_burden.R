# load packages

library(survival)
library(survminer)
library(dplyr)
library(readxl)
library(R.matlab)
library(openxlsx)

data_path<-"../../"
blca_data<-read_excel(paste(data_path,"blca_TMB.xlsx",sep=""),sheet="Sheet1")
tmb<-blca_data$`TMB values`
LI<-0.33
IH<-0.66
tLI<-quantile(tmb,probs=LI,names=FALSE)
tIH<-quantile(tmb,probs=IH,names=FALSE)
TMB_class<-cut(blca_data$`TMB values`,
               breaks=c(-1,tLI,tIH,Inf),labels=c("Low","Mid","High"))


blca_tcga<-read_excel(paste(data_path,"Table_S1.2017_08_05.xlsx",sep=""),sheet="Master table")

patient_ku<-blca_data$'Case ID'
patient_tcga<-blca_tcga$'Case ID'
tcga_futime<-blca_tcga$`Combined days to last followup or death`
tcga_fustat<-blca_tcga$'Vital status'

futime<-rep(0,length(patient_ku))
fustat<-character(length(patient_ku))
for (nn in 1:length(patient_ku)){
  pID_ss<-substring(patient_ku[nn],1,12)
  for (mm in 1:length(patient_tcga)){
    pID_oo<-as.vector(patient_tcga[mm])
    if (identical(pID_ss,pID_oo)){
      futime[nn]<-tcga_futime[mm]
      fustat[nn]<-tcga_fustat[mm]
      break
    }
  }
}


blca<-data.frame("patientID"=blca_data$`Case ID`,
                 "futime"=futime,
                 "fustat"=fustat,
                 "TMB_class"=TMB_class)

#write.xlsx(blca,"blca.xlsx")

# dichotomize the dead/alive and change data labels
# dead: censored 1, alive 0
blca$fustat <- factor(blca$fustat,
                      levels=c("Dead","Alive"),
                      labels=c("1","0"))


futime_d<-as.numeric(as.character(blca$futime))
row_ind<-which(is.na(futime_d))
blca<-blca[-c(row_ind),]   # remove two not available patients
blca$futime<-as.numeric(as.character(blca$futime))
blca$fustat<-as.numeric(as.character(blca$fustat))

# fit survival data using the kaplan-Meier method
surv_object<-Surv(time=blca$futime,event=blca$fustat)
# surv_object
# 
fit1<-survfit(surv_object~TMB_class,data=blca)
# summary(fit1)
# 
ggsurvplot(fit1,data=blca,pval=TRUE)

