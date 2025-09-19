# read in the exposure data
exposure.data <- read.csv(file = "Cmac-Vc-heat-inactivated-input-data.csv",stringsAsFactors = T)
names(exposure.data)
names(exposure.data)[names(exposure.data) == "cells.suspension"] <- "initial.damage"

# consider only surely NON-CONTAMINATED replicates and not partial replicate (without 10^6 CFU)
exposure.data <- exposure.data[exposure.data$Date.start != "2023-01-25",]
exposure.data <- exposure.data[exposure.data$Date.start != "2023-05-24",]

exposure.data <- droplevels(exposure.data)


#### Cox survival analyses ####
library(coxme)
library(survival)
library(ggplot2)
library(survminer)
library(emmeans)

### data prep: for each individual: age at death or last seen and status (censored or not)
# sum dead values for each individual
d.sum <- aggregate(exposure.data$Dead,by = list(id = exposure.data$Individual, Bacteria = exposure.data$Bacteria, CFU = exposure.data$CFU, Date.start = exposure.data$Date.start, initial.damage = exposure.data$initial.damage),sum)
# create a new dataframe with one row per individual
exposure.aad <- data.frame(d.sum[,c(1,2,3,4,5)])

# give censored status 0 to all the ones seen alive at 48h (sum of dead observations 0), 1 to all the others
exposure.aad$censored <-  1
exposure.aad[d.sum$x == 0,]$censored <- 0
# consider CFU as a categorical variable
exposure.aad$CFU <- as.factor(exposure.aad$CFU)

# add last observation (death or alive for the ones at 48h)
exposure.aad$last.obs <- NA
exposure.aad[d.sum$x == 0 | d.sum$x == 1,]$last.obs <- 48
exposure.aad[d.sum$x == 2,]$last.obs <- 24
exposure.aad[d.sum$x == 3,]$last.obs <- 6
exposure.aad[d.sum$x == 4,]$last.obs <- 2
exposure.aad[d.sum$x == 5,]$last.obs <- 1


hist(exposure.aad$last.obs)

### Cox mixed effects model
# null model
null.model <- survival::coxph(Surv(exposure.aad$last.obs,exposure.aad$censored)~1)
summary(null.model)


# basic Cox's proportional hazard model:
# time of last observation with censored status (survival object) as a function of CFU
# CFU is considered a discrete variable (factor), not continuous
base.model <- survival::coxph(Surv(exposure.aad$last.obs,exposure.aad$censored)~exposure.aad$CFU)
summary(base.model)

# Call:
#   survival::coxph(formula = Surv(exposure.aad$last.obs, exposure.aad$censored) ~ 
#                     exposure.aad$CFU)
# 
# n= 144, number of events= 16 
# 
# coef exp(coef) se(coef)      z Pr(>|z|)
# exposure.aad$CFU1000000 -0.2478    0.7806   0.5040 -0.492    0.623
# 
# exp(coef) exp(-coef) lower .95 upper .95
# exposure.aad$CFU1000000    0.7806      1.281    0.2907     2.096
# 
# Concordance= 0.53  (se = 0.064 )
# Likelihood ratio test= 0.24  on 1 df,   p=0.6
# Wald test            = 0.24  on 1 df,   p=0.6
# Score (logrank) test = 0.24  on 1 df,   p=0.6

car::Anova(base.model)

# Analysis of Deviance Table
# Cox model: response is Surv(exposure.aad$last.obs, exposure.aad$censored)
# Terms added sequentially (first to last)
# 
# loglik  Chisq Df Pr(>|Chi|)
# NULL              -78.652                     
# exposure.aad$CFU -78.530 0.2436  1     0.6216

# initial damage as explanatory variable
model.1 <- survival::coxph(Surv(exposure.aad$last.obs,exposure.aad$censored)~exposure.aad$CFU + exposure.aad$initial.damage)

summary(model.1)
# Call:
#   survival::coxph(formula = Surv(exposure.aad$last.obs, exposure.aad$censored) ~ 
#                     exposure.aad$CFU + exposure.aad$initial.damage)
# 
# n= 144, number of events= 16 
# 
# coef exp(coef) se(coef)      z Pr(>|z|)
# exposure.aad$CFU1000000     -0.3223    0.7245   0.5133 -0.628    0.530
# exposure.aad$initial.damage  1.4501    4.2636   1.0526  1.378    0.168
# 
# exp(coef) exp(-coef) lower .95 upper .95
# exposure.aad$CFU1000000        0.7245     1.3802    0.2649     1.981
# exposure.aad$initial.damage    4.2636     0.2345    0.5417    33.556
# 
# Concordance= 0.535  (se = 0.069 )
# Likelihood ratio test= 1.55  on 2 df,   p=0.5
# Wald test            = 2.04  on 2 df,   p=0.4
# Score (logrank) test = 2.3  on 2 df,   p=0.3

car::Anova(model.1)
# Analysis of Deviance Table (Type II tests)
# LR Chisq Df Pr(>Chisq)
# exposure.aad$CFU             0.39865  1     0.5278
# exposure.aad$initial.damage  1.30878  1     0.2526

anova(null.model,model.1)
# Analysis of Deviance Table
# Cox model: response is  Surv(exposure.aad$last.obs, exposure.aad$censored)
# Model 1: ~ 1
# Model 2: ~ exposure.aad$CFU + exposure.aad$initial.damage
# loglik  Chisq Df Pr(>|Chi|)
# 1 -78.652                     
# 2 -77.876 1.5524  2     0.4602


######################
### Plots ###########
######################
# plot base model
ggsurv <- survminer::ggsurvplot(survfit(Surv(last.obs,censored)~CFU,data=exposure.aad),data=exposure.aad,
                                conf.int = TRUE,censor = FALSE, palette="grey",
                                legend.title = "Bacterial Load", legend.labs = c("Control", "High"))
ggsurv$plot + theme_classic()+
  theme(panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))+
  scale_y_continuous(limits = c(0.5,1))+
  scale_x_continuous(breaks = seq(0,50,5))+
  labs(x = "Exposure hours")

ggsave("grey_Cmac_adults_Vc_heat-inactivated_Vc_survplot.pdf", device = "pdf", path = ".",
       dpi = 300, width = 150, height = 105, units = "mm")
