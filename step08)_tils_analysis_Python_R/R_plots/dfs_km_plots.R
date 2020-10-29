# dfs survival plots for the paper

library(survival)
library(survminer)
library(dplyr)
library(readxl)
library(openxlsx)

dfs<-read_excel(paste("./","TCGA_BLCA_Fig.5_dfs.xlsx",sep=""),sheet="Sheet1")

figa<-TRUE
figb<-TRUE
figc<-TRUE
figd<-TRUE

# row_ind<-which(is.na(dfs$`Disease Free (Months)`))
# if (length(row_ind)>0) {
#   dfs<-dfs[-c(row_ind),]   # remove two not available patients
# }

# row_ind2<-which(dfs$`Disease Free (Months)`<0)
# if (length(row_ind2)>0){
#   dfs<-dfs[-c(row_ind2),]
# }
# 
# row_ind3<-which(dfs$`Disease Free Status`=='NA')
# if (length(row_ind3)>0){
#   dfs<-dfs[-c(row_ind3),]
# }


if (figa==TRUE){
  row_ind4<-which(dfs$figa=='NA')
  if (length(row_ind4)>0){
    dfs2<-dfs[-c(row_ind4),]
  }
  
  dfs2$fustat <- factor(dfs2$`Disease Free Status`,
                       levels=c("1:Recurred/Progressed","0:DiseaseFree"),
                       labels=c("1","0"))
  
  dfs2$futime<-as.numeric(as.character(dfs2$`Disease Free (Months)`)) # note that: must be numeric type
  dfs2$fustat<-as.numeric(as.character(dfs2$fustat)) # note that: must be numeric type
  
  surv_object<-Surv(time=dfs2$futime,event=dfs2$fustat)
  
  fit1<-survfit(surv_object~figa,data=dfs2)
  
  survp<-ggsurvplot(fit1,pval = TRUE,
             risk.table = TRUE, 
             legend=c(0.75,0.85),
             #legend.labs=c("High-High-High-Low (HHHL) (27)", "Low-High-High-Low (LHHL) (11)","High (w/o HHHL, LHHL) (88)"),
             legend.title="Categories",
             xlab="Time in months")+ggtitle("TCGA BLCA")
  
  ggsave(file = "./figa_dfs.pdf", print(survp,newpage = FALSE))
}

if (figb==TRUE){
  row_ind4<-which(dfs$figb=='NA')
  if (length(row_ind4)>0){
    dfs2<-dfs[-c(row_ind4),]
  }
  
  dfs2$fustat <- factor(dfs2$`Disease Free Status`,
                       levels=c("1:Recurred/Progressed","0:DiseaseFree"),
                       labels=c("1","0"))
  
  dfs2$futime<-as.numeric(as.character(dfs2$`Disease Free (Months)`)) # note that: must be numeric type
  dfs2$fustat<-as.numeric(as.character(dfs2$fustat)) # note that: must be numeric type
  
  surv_object<-Surv(time=dfs2$futime,event=dfs2$fustat)
  
  fit1<-survfit(surv_object~figb,data=dfs2)
  
  survp<-ggsurvplot(fit1,pval = TRUE,
                    risk.table = TRUE, 
                    legend=c(0.75,0.85),
                    #legend.labs=c("High-High-High-Low (HHHL) (27)", "Low-High-High-Low (LHHL) (11)","High (w/o HHHL, LHHL) (88)"),
                    legend.title="Categories",
                    xlab="Time in months")+ggtitle("TCGA BLCA")
  
  ggsave(file = "./figb_dfs.pdf", print(survp,newpage = FALSE))
}

if (figc==TRUE){
  row_ind4<-which(dfs$figc=='NA')
  if (length(row_ind4)>0){
    dfs2<-dfs[-c(row_ind4),]
  }
  
  dfs2$fustat <- factor(dfs2$`Disease Free Status`,
                       levels=c("1:Recurred/Progressed","0:DiseaseFree"),
                       labels=c("1","0"))
  
  dfs2$futime<-as.numeric(as.character(dfs2$`Disease Free (Months)`)) # note that: must be numeric type
  dfs2$fustat<-as.numeric(as.character(dfs2$fustat)) # note that: must be numeric type
  
  surv_object<-Surv(time=dfs2$futime,event=dfs2$fustat)
  
  fit1<-survfit(surv_object~figc,data=dfs2)
  
  survp<-ggsurvplot(fit1,pval = TRUE,
                    risk.table = TRUE, 
                    legend=c(0.75,0.85),
                    #legend.labs=c("High-High-High-Low (HHHL) (27)", "Low-High-High-Low (LHHL) (11)","High (w/o HHHL, LHHL) (88)"),
                    legend.title="Categories",
                    xlab="Time in months")+ggtitle("TCGA BLCA")
  
  ggsave(file = "./figc_dfs.pdf", print(survp,newpage = FALSE))
}

if (figd==TRUE){
  row_ind4<-which(dfs$figd=='NA')
  if (length(row_ind4)>0){
    dfs2<-dfs[-c(row_ind4),]
  }
  
  dfs2$fustat <- factor(dfs2$`Disease Free Status`,
                        levels=c("1:Recurred/Progressed","0:DiseaseFree"),
                        labels=c("1","0"))
  
  dfs2$futime<-as.numeric(as.character(dfs2$`Disease Free (Months)`)) # note that: must be numeric type
  dfs2$fustat<-as.numeric(as.character(dfs2$fustat)) # note that: must be numeric type
  
  surv_object<-Surv(time=dfs2$futime,event=dfs2$fustat)
  
  fit1<-survfit(surv_object~figd,data=dfs2)
  
  survp<-ggsurvplot(fit1,pval = TRUE,
                    risk.table = TRUE, 
                    legend=c(0.75,0.85),
                    #legend.labs=c("High-High-High-Low (HHHL) (27)", "Low-High-High-Low (LHHL) (11)","High (w/o HHHL, LHHL) (88)"),
                    legend.title="Categories",
                    xlab="Time in months")+ggtitle("TCGA BLCA")
  
  ggsave(file = "./figd_dfs.pdf", print(survp,newpage = FALSE))
}
