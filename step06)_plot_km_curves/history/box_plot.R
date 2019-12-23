library(R.matlab)
library(readxl)

#ddf = data.frame(NUMS = rnorm(500), GRP = sample(LETTERS[1:5],500,replace=T))
#boxplot(NUMS ~ GRP, data = ddf, lwd = 2, ylab = 'NUMS')
#stripchart(NUMS ~ GRP, vertical = TRUE, data = ddf, 
#           method = "jitter", add = TRUE, pch = 20, col = 'blue')

data_path0<-"E:/Hongming/projects/tcga-bladder-mutationburden/heatmap_blca/"
pid_pred_gt<-read_excel(paste(data_path0,"pid_pred_gt.xlsx",sep=""),sheet = "Sheet1")

path <- system.file("mat-files", package = "R.matlab")
pathname <- file.path("E:/Hongming/projects/tcga-bladder-mutationburden/Hongming_codes/step11)_tmb_map_process/features/", "feat_hl.mat")
blca_pre<-readMat(pathname)

pid<-vector()
plabel<-vector()
pvalue_h<-vector()
pvalue_l<-vector()
pvalue<-vector()

for (kk in 1:length(blca_pre[[1]])/3){
  pid[kk]<-blca_pre[[1]][[kk]]
  pvalue_h[kk]<-as.numeric(blca_pre[[1]][[kk+length(blca_pre[[1]])/3]])
  pvalue_l[kk]<-as.numeric(blca_pre[[1]][[kk+length(blca_pre[[1]])*2/3]])
  
}

tt=0
for (kk in 1:length(pid)){
  temp_pID<-substring(as.character(pid[kk]),1,23)
  for (jj in 1:length(pid_pred_gt$`Patient ID`))
    if (temp_pID==substring(as.character(pid_pred_gt$`Patient ID`[jj]),2,24)){
      if (pid_pred_gt$Ground_Truths[jj]=="'High'" & pid_pred_gt$Automatic_Predictions[jj]=="'High'"){
        plabel[kk]<-paste('High','-','True',sep="")
        tt=tt+1
        pvalue[kk]<-pvalue_h[kk]/pvalue_l[kk]
        
      } else if (pid_pred_gt$Ground_Truths[jj]=="'Low'" & pid_pred_gt$Automatic_Predictions[jj]=="'Low'"){
        plabel[kk]<-paste('Low','-','True',sep="")
        pvalue[kk]<-pvalue_l[kk]/pvalue_h[kk]
        tt=tt+1
        
      } else if (pid_pred_gt$Ground_Truths[jj]=="'Low'" & pid_pred_gt$Automatic_Predictions[jj]=="'High'"){
        plabel[kk]<-paste('High','-','False',sep="")
        pvalue[kk]<-pvalue_h[kk]/pvalue_l[kk]
        
      } else if (pid_pred_gt$Ground_Truths[jj]=="'High'" & pid_pred_gt$Automatic_Predictions[jj]=="'Low'"){
        plabel[kk]<-paste('Low','-','False',sep="")
        pvalue[kk]<-pvalue_l[kk]/pvalue_h[kk]
        
      } else {
        print('not possible!!')
      }
          
      break
    }
}

#temp=as.numeric(pvalue)
#max_value=max(temp[temp!=Inf])
pvalue[pvalue>20]<-20


ddf=data.frame(proportion=as.numeric(pvalue),label=plabel)

#ddf = data.frame(NUMS = rnorm(500), GRP = sample(LETTERS[1:5],500,replace=T))
boxplot(proportion ~ label, data = ddf, lwd = 2, ylab = 'proportion')
stripchart(proportion ~ label, vertical = TRUE, data = ddf, 
           method = "jitter", add = TRUE, pch = 20, col = 'blue')

#ggplot(ddf, aes(x=plabel, y=proportion)) + 
#  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
#  geom_jitter(position=position_jitter(width=.1, height=0))