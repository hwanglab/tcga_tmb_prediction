# good reference website:http://www.sthda.com/english/wiki/cox-proportional-hazards-model

library("survminer");
library(survival)

library(ggplot2)

#setwd("D:\\Research\\Code\\fromHongming\\")

#--- data preprocessing 
# load the tabel 
mstr_file_name = paste(".\\",
                       "Table_S1.2017_08_05_predictions.txt", sep="", collapse = "");
my_data <- read.table(mstr_file_name, sep="\t", header=TRUE);
my_data$Combined.days.to.last.followup.or.death = as.numeric(my_data$Combined.days.to.last.followup.or.death)/30

old_vitalstatus = my_data$Vital.status

my_data$Vital.status[old_vitalstatus=="Dead"] = 1
my_data$Vital.status[old_vitalstatus=="Alive"] = 0

mv_selIDX = !is.na(my_data$Combined.days.to.last.followup.or.death)
my_data_sub = my_data[mv_selIDX, ]

my_data_sub$Vital.status = as.numeric(my_data_sub$Vital.status)
my_data_sub[my_data_sub=="ND"]<-NA

# use "factor" for the categorical variables  
my_data_sub$AJCC.pathologic.tumor.stage = factor(my_data_sub$AJCC.pathologic.tumor.stage, levels = c( "II", "III", "IV" ))
my_data_sub$Inflammatory.Infiltrate.Response = factor(my_data_sub$Inflammatory.Infiltrate.Response)
my_data_sub$Lymphovascular.invasion = factor(my_data_sub$Lymphovascular.invasion)

my_data_sub$RPPA.cluster = factor(my_data_sub$RPPA.cluster)
my_data_sub$mRNA.cluster = factor(my_data_sub$mRNA.cluster)

my_data_sub$fb = factor(my_data_sub$fb)


#--- Multivariate cox model
fit <- coxph(formula = Surv(Combined.days.to.last.followup.or.death, Vital.status) 
             ~ Age.at.diagnosis + AJCC.pathologic.tumor.stage + Inflammatory.Infiltrate.Response 
             + Lymphovascular.invasion + RPPA.cluster + mRNA.cluster + fb, data = my_data_sub)

#summary(fit)

# save texts (on the command windows) to a text file   
sink(file='./cox_multivariate_summary.txt')

print(summary(fit))  #output provides HR CIs
print(confint(fit))  #coefficient CIs
print(exp(confint(fit)))  #Also HR CIs

sink() # set the output back to the console
#sink(NULL)

# Prepare the columns
beta <- exp(coef(fit))  ## hazard ratio
CI   <- round(exp(confint(fit)), 5) ## hazard ratio of confidence intervals

coeffs <- coef(summary(fit))
p    <- as.matrix(coeffs[,5])   ## the 5 th column is the p-values

# Bind columns together, and select desired rows
res <- cbind(beta, CI, p)

fileConn<-file('./cox_multivariate_CIs.txt', "w")
for (mn_i in 1:length(beta)){
  mstr_line = sprintf("%s\t %.5f (%.5f, %.5f)\t %.5f", rownames(res)[mn_i],
                      beta[mn_i], CI[mn_i,1], CI[mn_i,2], p[mn_i])
  write(mstr_line, fileConn, append=T) 
}
close(fileConn)  





