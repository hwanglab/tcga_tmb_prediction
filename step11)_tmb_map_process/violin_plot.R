# load packages

library(survival)
library(survminer)
library(dplyr)
library(readxl)
library(R.matlab)
library(openxlsx)
library(ggplot2)


output_path<-"E:/Hongming/projects/tcga-bladder-mutationburden/heterogeneity/violin_log_1/"
#output_path<-'E:/Hongming/projects/tcga-bladder-mutationburden/heterogeneity/violin_zscore_1/'

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
      plabel[kk]<-paste(pid_pred_gt$Automatic_Predictions[jj],'-',plabel1[kk],sep="")
      break
    }
}

pid2<-vector()
plabel2<-vector()
ind=1
for (kk in 1:length(pid)){
  if (plabel[kk]=="'High'-Low") {
    pid2[ind]<-pid[kk]
    plabel2[ind]<-'High-Low'
    ind<-ind+1
  } else if (plabel[kk]=="'High'-High" ){
    pid2[ind]<-pid[kk]
    plabel2[ind]<-'High-High'
    ind<-ind+1
  } else if (plabel[kk]=="'Low'-High"){
    pid2[ind]<-pid[kk]
    plabel2[ind]<-'Low-High'
    ind<-ind+1
  } else if (plabel[kk]=="'Low'-Low"){
    pid2[ind]<-pid[kk]
    plabel2[ind]<-'Low-Low'
    ind<-ind+1
  } 
  else {
    print('skip!')
  }
  
}

data_path1<-"E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/step11)_tmb_map_process/"
pid_pdl1<-read_excel(paste(data_path1,"cBioPortal_data-3.xlsx",sep=""),sheet = "cBioPortal_data-3")


for (cc in 1:length(pid_pdl1$COMMON)){
#for (cc in 13:13){                      # consider only cd274
#cc=13
  ppdl1<-vector()
  for (kk in 1:length(pid2)){
    pid_temp<-as.character(pid2[kk])
    pid_samp<-paste(substring(pid_temp,1,12),'-01',sep="")
    if (!is.null(pid_pdl1[[pid_samp]][cc])){
      ppdl1[kk]<-log(as.numeric(pid_pdl1[[pid_samp]][cc])) # CD274
    }
  }
  
  ddf=data.frame(PDL1=as.numeric(ppdl1),label=plabel2)
  #boxplot(PDL1 ~ label, data = ddf, lwd = 2, ylab = 'pdl1')
  #stripchart(proportion ~ label, vertical = TRUE, data = ddf, 
  #           method = "jitter", add = TRUE, pch = 20, col = 'blue')
  
  
  
  # convert the variable dose from a numeric to a factor varibale
  #ToothGrowth$dose<-as.factor(ToothGrowth$dose)
  #head(ToothGrowth)
  
  #setEPS()
  #postscript("whatever.eps")
  
  # basic violin plot
  p<-ggplot(ddf,aes(x=label,y=PDL1))+geom_violin(trim=FALSE,fill="grey90",color="darkred")+
    labs(title="Violin plot of gene expression",x="categories", y = pid_pdl1$COMMON[cc])+
    geom_boxplot(width=0.1)+
    geom_jitter(shape=16, position=position_jitter(0.2))+ stat_summary(fun.data=mean_sdl, fun.args=list(mult=1), 
                                    geom="pointrange", color="red")
  
  p+stat_compare_means(method="aov")
  
  #dev.off()
  
  # violin plot with mean points
  #p + stat_summary(fun.y=mean, geom="point", shape=23, size=2)
  # violin plot with median points
  #p + stat_summary(fun.y=median, geom="point", size=2, color="red")
  #p + geom_boxplot(width=0.1)
  #p + stat_summary(fun.data="mean_sdl", mult=1, 
  #                 geom="crossbar", width=0.2 )
  
  ggsave(paste(output_path,pid_pdl1$COMMON[cc],'.png'))
}
