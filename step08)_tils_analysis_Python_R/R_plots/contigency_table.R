# make contigency table

library(survival)
library(survminer)
library(dplyr)
library(readxl)


data_path0<-"E:/Hongming/papers/2020 cancer patient tmb prediction/TMB_Manuscript/"
df<-read_excel(paste(data_path0,"TCGA_BLCA_Fig.5.xlsx",sep=""))


plabel0<-vector()
pvalue<-df$'Tils_density'
tt<-quantile(as.numeric(unlist(pvalue)),c(0.333,0.5,0.667),na.rm = TRUE)
plabel0<-(pvalue>tt[2])
plabel0[plabel0==TRUE]<-'High'
plabel0[plabel0==FALSE]<-'Low'

plabel1<-df$'Pred_tmb'

plabel2<-vector()
pvalue<-df$'Pred_entropy'
tt<-quantile(as.numeric(unlist(pvalue)),c(0.333,0.5,0.667),na.rm = TRUE)
plabel2<-(pvalue>tt[2])
plabel2[plabel2==TRUE]<-'High'
plabel2[plabel2==FALSE]<-'Low'
plabel2[which(pvalue=='NA')]<-'NA'

row_ind<-which(is.na(df$'Combined days to last followup or death'))
if (length(row_ind)>0) {
  plabel2[row_ind]<-'NA'
}

row_ind2<-which(df$'Combined days to last followup or death'<0)
if (length(row_ind)>0){
  plabel2[row_ind2]<-'NA'
}

hhl<-0
hhh<-0
hll<-0
hlh<-0

lhl<-0
lhh<-0
lll<-0
llh<-0
tt<-0

for (i in 1:length(plabel0))
{
  if (plabel0[i]=='High' && plabel1[i]=='High' && plabel2[i]=='Low'){
    hhl<-hhl+1
  } else if (plabel0[i]=='High' && plabel1[i]=='High' && plabel2[i]=='High'){
    hhh<-hhh+1
  } else if (plabel0[i]=='High' && plabel1[i]=='Low' && plabel2[i]=='Low'){
    hll<-hll+1
  } else if (plabel0[i]=='High' && plabel1[i]=='Low' && plabel2[i]=='High'){
    hlh<-hlh+1
  } else if (plabel0[i]=='Low' && plabel1[i]=='High' && plabel2[i]=='Low'){
    lhl<-lhl+1
  } else if (plabel0[i]=='Low' && plabel1[i]=='High' && plabel2[i]=='High'){
    lhh<-lhh+1
  } else if (plabel0[i]=='Low' && plabel1[i]=='Low' && plabel2[i]=='Low'){
    lll<-lll+1
  } else if (plabel0[i]=='Low' && plabel1[i]=='Low' && plabel2[i]=='High'){
    llh<-llh+1
  } else{
    print(i)
    tt<-tt+1
  }
}

t=0
