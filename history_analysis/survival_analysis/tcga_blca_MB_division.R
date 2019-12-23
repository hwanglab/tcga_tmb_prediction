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
tmb<-tmb/50
tLI<-5

TMB_class<-cut(tmb,breaks=c(-1,tLI,Inf),labels=c("Low","Midhig"))

blca<-data.frame("futime"=blca_data$`Combined days to last followup or death`,
                 "fustat"=blca_data$`Vital status`,
                 "TMB"=blca_data$`total number Single Nucleotide Variants`,
                 "TMB_class"=TMB_class)

#write.xlsx(blca,"blca.xlsx")

# dichotomize the dead/alive and change data labels
# dead: censored 1, alive 0
blca$fustat <- factor(blca$fustat,
                      levels=c("Dead","Alive"),
                      labels=c("1","0"))

futime_d<-as.numeric(as.character(blca$futime))
row_ind<-which(is.na(futime_d))
blca2<-blca[-c(row_ind),]   # remove two not available patients
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