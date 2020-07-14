library("survival")
library("survminer")
data("lung")
head(lung)
res.cox <- coxph(Surv(time, status) ~ pat.karno, data = lung)
res.cox
