

# onw-way ANOVA test
# ref: http://www.sthda.com/english/wiki/one-way-anova-test-in-r

my_data<-PlantGrowth

# shown a random sample
set.seed(1234)
dplyr::sample_n(my_data,10)

# show the levels
levels(my_data$group)

library(dplyr)
group_by(my_data,group) %>%
  summarise(
    count=n(),
    mean=mean(weight,na.rm = TRUE),
    sd=sd(weight,na.rm = TRUE)
  )

# Box plots
# ++++++++++++++++++++
# Plot weight by group and color by group
library("ggpubr")
ggboxplot(my_data, x = "group", y = "weight", 
          color = "group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("ctrl", "trt1", "trt2"),
          ylab = "Weight", xlab = "Treatment")


# Mean plots
# ++++++++++++++++++++
# Plot weight by group
# Add error bars: mean_se
# (other values include: mean_sd, mean_ci, median_iqr, ....)
library("ggpubr")
ggline(my_data, x = "group", y = "weight", 
       add = c("mean_se", "jitter"), 
       order = c("ctrl", "trt1", "trt2"),
       ylab = "Weight", xlab = "Treatment")


# compute the analysis of variance
res.aov<-aov(weight~group,data=my_data)
# summary of the analysis
summary(res.aov)