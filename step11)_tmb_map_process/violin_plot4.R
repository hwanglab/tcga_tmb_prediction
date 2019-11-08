# load packages

library(survival)
library(survminer)
library(dplyr)
library(readxl)
library(R.matlab)
library(openxlsx)
library(ggplot2)


output_path<-"E:/Hongming/projects/tcga-bladder-mutationburden/heterogeneity/gene_plots/wex_log/"
excel_output<-"E:/Hongming/projects/tcga-bladder-mutationburden/heterogeneity/gene_plots/"

data_path0<-"E:/Hongming/projects/tcga-bladder-mutationburden/heatmap_blca/"
pid_pred_gt<-read_excel(paste(data_path0,"pid_pred_gt.xlsx",sep=""),sheet = "Sheet1")
pid<-pid_pred_gt$`Patient ID`

plabel<-vector()
for (kk in 1:length(pid)){
  temp_pID<-substring(as.character(pid[kk]),2,24)
  plabel[kk]<-pid_pred_gt$Ground_Truths[kk]
}

pid2<-vector()
plabel2<-vector()
ind=1
for (kk in 1:length(pid)){
  if (plabel[kk]=="'High'") {
    pid2[ind]<-substring(as.character(pid[kk]),2,24)
    plabel2[ind]<-'High'
    ind<-ind+1
  } else if (plabel[kk]=="'Low'"){
    pid2[ind]<-substring(as.character(pid[kk]),2,24)
    plabel2[ind]<-'Low'
    ind<-ind+1
  } 
  else {
    print('skip!')
  }
  
}

data_path1<-"E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/step11)_tmb_map_process/"
pid_pdl1<-read_excel(paste(data_path1,"cBioPortal_data-3.xlsx",sep=""),sheet = "cBioPortal_data-3")

#xx=data.frame(GENE_ID=pid_pdl1$GENE_ID,COMMON=pid_pdl1$COMMON)
gene_id<-vector()
common<-vector()
pvalue<-vector()

for (cc in 1:length(pid_pdl1$COMMON)){
  #for (cc in 13:13){                      # consider only cd274
  #cc=13
  gene_id[cc]<-pid_pdl1$GENE_ID[cc]
  common[cc]<-pid_pdl1$COMMON[cc]
  
  ppdl1<-vector()
  for (kk in 1:length(pid2)){
    pid_temp<-as.character(pid2[kk])
    pid_samp<-paste(substring(pid_temp,1,12),'-01',sep="")
    if (!is.null(pid_pdl1[[pid_samp]][cc])){
      ppdl1[kk]<-log10(as.numeric(pid_pdl1[[pid_samp]][cc])+10) # CD274
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
  #my_comparisons<-list(c('High-Low',"Low-Low"))
  # basic violin plot
  p<-ggplot(ddf,aes(x=label,y=PDL1))+geom_violin(trim=FALSE,fill="grey90",color="darkred")+
    labs(title="Violin plot of gene expression",x="categories", y = pid_pdl1$COMMON[cc])+
    geom_boxplot(width=0.1)+
    geom_jitter(shape=16, position=position_jitter(0.2))+ 
    stat_summary(fun.data=mean_sdl, fun.args=list(mult=1), geom="pointrange", color="red")
    #scale_y_continuous(trans='log10')
    
  
  p+stat_compare_means(method="aov") 
  #stat_compare_means(label.y = 50)     # Add global p-value, at position y=50
  ptemp<-compare_means(PDL1~label,ddf,method="aov")
  
  pvalue[cc]<-ptemp$p
  
  #dev.off()
  
  # violin plot with mean points
  #p + stat_summary(fun.y=mean, geom="point", shape=23, size=2)
  # violin plot with median points
  #p + stat_summary(fun.y=median, geom="point", size=2, color="red")
  #p + geom_boxplot(width=0.1)
  #p + stat_summary(fun.data="mean_sdl", mult=1, 
  #                 geom="crossbar", width=0.2 )
  
  ggsave(paste(output_path,pid_pdl1$COMMON[cc],sprintf("_pvalue=%f_",ptemp$p),'.png',sep = ""))
}

## write excel file
xx=data.frame(GENE_ID=gene_id,COMMON=common, WEX_High_Low_pvalue=pvalue)
write.xlsx(xx,paste(excel_output,'gene_log10_pvalues.xlsx'),sheetName="Sheet1",
           col.names=TRUE,append=FALSE)