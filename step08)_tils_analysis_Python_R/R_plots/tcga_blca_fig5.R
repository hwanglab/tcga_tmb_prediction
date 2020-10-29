#generate the master table for fig.5 in tcga blca

library(survival)
library(survminer)
library(dplyr)
library(readxl)
library(openxlsx)

dfs<-read_excel(paste("./","2017_blca_clinical_info.xlsx",sep=""),sheet="Sheet1")

fis<-read_excel(paste("./","TCGA_BLCA_Fig.5.xlsx",sep = ""),sheet="Sheet1")


dfs_time<-vector()
dfs_status<-vector()
os_time<-vector()
os_status<-vector()
for (nn in 1 :length(fis$`Case ID`)){
   temp_pID=as.character(fis$`Case ID`[nn])
   ind<-which(dfs$`Patient ID`==temp_pID)
   if (length(ind)==1) {
     dfs_time<-c(dfs_time,as.numeric(dfs$`Disease Free (Months)`[ind]))
     dfs_status<-c(dfs_status,as.character(dfs$`Disease Free Status`[ind]))
     
     os_time<-c(os_time,as.numeric(dfs$'Overall Survival (Months)'[ind]))
     os_status<-c(os_status,as.character(dfs$'Overall Survival Status'[ind]))
              
   } else {
     dfs_time<-c(dfs_time,'NA')
     dfs_status<-c(dfs_status,'NA')
     
     os_time<-c(os_time,'NA')
     os_status<-c(os_status,'NA')
   }
 }
 
fis['Disease Free (Months)']<-dfs_time
fis['Disease Free Status']<-dfs_status

fis['Overall Survival (Months)']<-os_time
fis['Overall Survival Status']<-os_status

write.xlsx(fis, paste('./',"TCGA_BLCA_Fig.5_dfs.xlsx",sep=""), sheetName = "Sheet1")