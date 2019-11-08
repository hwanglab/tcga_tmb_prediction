
library("readxl")
library(dplyr)
library("ggpubr")

pid<-read_excel("E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/pathway/TCGA_Data_BLCA_10042019.xlsx",sheet = 'PatientIDs',col_names = FALSE)
RNA_genes<-read_excel("E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/pathway/TCGA_Data_BLCA_10042019.xlsx",sheet = 'RNA_genes',col_names = FALSE)
RNA_X<-read.csv(file="E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/pathway/TCGA_Data_BLCA_10042019.csv",sep=',',header = FALSE)


excel_output<-"E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/pathway/pathway_plots/"
pathway<-read_excel("E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/pathway/kegg_pathway.xlsx",sheet='OrgDocu',col_names = FALSE)

data_path0<-"E:/Hongming/projects/tcga-bladder-mutationburden/heatmap_blca/"
pid_pred_gt<-read_excel(paste(data_path0,"pid_pred_gt.xlsx",sep=""),sheet = "Sheet1")

pathway_name<-vector()
wex_pvalues<-vector()

for (i in 1: nrow(pathway)) # number of pathways
{
  x<-data.frame("patient ID"=pid[[1]])
  
  pathway_name[i]<-pathway$X__1[i]
  
  for (j in 3:ncol(pathway)) # membership gene on each pathway
  {
    ms_gene<-pathway[i,j]
    ind<-pmatch(ms_gene[[1]],RNA_genes[[1]])
    if (!is.na(ind)) {
      if (!any(is.na((RNA_X[ind]))))
      {
        x[ms_gene[[1]]]<-RNA_X[ind]
      }
    }
  }
  mev<-apply(x[c(2:ncol(x))], 1, mean)
  x['mev']<-log10(mev+10)
  
  plabel<-rep('unknown',nrow(x))
  for (k in 1:nrow(x))
  {
    pid_temp<-substr(x[[1]][k],2,13) # patient id
    
    for (jj in 1:length(pid_pred_gt$`Patient ID`))
      if (pid_temp==substring(as.character(pid_pred_gt$`Patient ID`[jj]),2,13)){
        temp<-pid_pred_gt$Ground_Truths[jj]
        plabel[k]<-substr(temp,2,nchar(temp)-1)
        break
      }
  }
  x['plabel']<-factor(plabel)
  x2<-x[c(1,ncol(x)-1,ncol(x))]
  row_ind<-which(x2$plabel=='unknown')
  x3<-x2[-c(row_ind),]
  
  group_by(x3,plabel) %>%
    summarise(
      count=n(),
      mean=mean(mev,na.rm=TRUE),
      sd=sd(mev,na.rm = TRUE)
    )
  
 #p<-ggboxplot(x3,x="plabel",y="mev",
#            color="plabel",palette = c("#00AFBB", "#E7B800"),
#            order=c("High","Low"),
#           ylab="mean express values",xlab="predicted TMB groups")
  
# violin plot with dots
  p<-ggplot(x3,aes(x=plabel,y=mev))+geom_violin(trim = FALSE,fill="grey90",color="darkred")+
    labs(title="Plot of gene express values by tmb levels",x="WEX TMB groups", y = "mean express values")+geom_boxplot(width=0.1)+
    geom_jitter(shape=16, position=position_jitter(0.2))    
    
  
  p+stat_compare_means(method="aov")
  
  ptemp<-compare_means(mev~plabel,x3,method="aov")
  
  wex_pvalues[i]<-ptemp$p
  
  ggsave(paste('E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/pathway/pathway_plots/wex_pathway_plots/',
               pathway_name[i],sprintf("_pvalue=%f_",ptemp$p),'.png',sep=""),width = 10, height = 10)
  #p+stat_compare_means(method="t.test")
  
  # compute one-way anova test
  #res.aov<-aov(mev~plabel,data=x3)
  #pvalue<-summary(res.aov)[[1]][["Pr(>F)"]][[1]]
  #summary(res.aov)
  
}

## write excel file
xx=data.frame(Pathway_Name=pathway_name, WEX_High_Low_pvalue=wex_pvalues)
write.xlsx(xx,paste(excel_output,'pathway_log10_pvalues.xlsx',sep=""),sheetName="Sheet1",
           col.names=TRUE,append=FALSE)
