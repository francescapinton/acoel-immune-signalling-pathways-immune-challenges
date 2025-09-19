# read in the exposure data
exposure.data.init <- read.csv(file = "Cmac-Priestia-cleaned-up-data.csv",stringsAsFactors = T)
names(exposure.data.init)
names(exposure.data.init)[names(exposure.data.init) == "cells.suspension"] <- "initial.damage"

# removed first pilot replicate (very low numbers of 10^6 and control)
exposure.data <- exposure.data.init[exposure.data.init$Date.start != "2024-05-27",]
exposure.data <- droplevels(exposure.data)


###############################
#### Cox survival analyses ####
###############################

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
exposure.aad[d.sum$x == 4,]$last.obs <- 4
exposure.aad[d.sum$x == 5,]$last.obs <- 2
exposure.aad[d.sum$x == 6,]$last.obs <- 1

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

#### build and test models
anova(null.model,base.model)
# base model significantly better than null model

# Fixed effects: CFU, cell.suspension
# Random effects: experiment/replicate (Date.start), Individual

me.model <- coxme::coxme(Surv(last.obs,censored)~CFU*initial.damage + (1 | Date.start / id ),data = exposure.aad)
summary(me.model)

### simplify the model
# remove the (non-significant) interaction between CFU and initial.damage
me.model1 <- coxme::coxme(Surv(last.obs,censored)~CFU + initial.damage + (1 | Date.start / id ),data = exposure.aad)
summary(me.model1)


anova(me.model,me.model1)
# accept the simpler model me.model1

# no id as nested random effect
me.model2 <- coxme::coxme(Surv(last.obs,censored)~CFU + initial.damage + (1 | Date.start),data = exposure.aad)
summary(me.model2)

anova(me.model1,me.model2)
# p > 0.8 -> accept simpler model without id

# removed cells in suspension factor
me.model3 <- coxme::coxme(Surv(last.obs,censored)~CFU + (1 | Date.start ),data = exposure.aad)
summary(me.model3)
anova(me.model2,me.model3)
# p = 0.006 > 0.001 -> remove initial damage 

# remove random effects completely
fe.model <- survival::coxph(Surv(last.obs,censored)~CFU ,data = exposure.aad)
summary(fe.model)

anova(me.model2,fe.model)
# p << 0.001 -> random effect of experiment Start date is important

# compare base model to me.model3
anova(base.model,me.model3)
# me.model3 significantly better --> kept

# minimal adequate model is me.model3 coxme::coxme(Surv(last.obs,censored)~CFU + (1 | Date.start ),data = exposure.aad)
summary(me.model3)

# Cox mixed-effects model fit by maximum likelihood
# Data: exposure.aad
# events, n = 354, 644
# Iterations= 14 74 
# NULL Integrated    Fitted
# Log-likelihood -2167.331  -2032.316 -2027.967
# 
# Chisq   df p    AIC    BIC
# Integrated loglik 270.03 3.00 0 264.03 252.42
# Penalized loglik 278.73 3.93 0 270.87 255.66
# 
# Model:  Surv(last.obs, censored) ~ CFU + (1 | Date.start) 
# Fixed coefficients
# coef exp(coef)  se(coef)     z p
# CFU1e+05 -0.000490705 0.9995094 0.1601440  0.00 1
# CFU1e+06  1.859511804 6.4206015 0.1416685 13.13 0
# 
# Random effects
# Group      Variable  Std Dev   Variance 
# Date.start Intercept 0.5004184 0.2504186

## post hoc comparisons

me.model3.posthoc <-  emmeans(me.model3, specs = pairwise ~ CFU, type = "response",adjust = "bonferroni")
summary(me.model3.posthoc,infer = TRUE)

# $emmeans
# CFU   response     SE  df asymp.LCL asymp.UCL null z.ratio p.value
# 0        0.538 0.0478 Inf     0.452      0.64    1  -6.981  <.0001
# 1e+05    0.537 0.0480 Inf     0.451      0.64    1  -6.955  <.0001
# 1e+06    3.452 0.2692 Inf     2.962      4.02    1  15.884  <.0001
# 
# Confidence level used: 0.95 
# Intervals are back-transformed from the log scale 
# Tests are performed on the log scale 
# 
# $contrasts
# contrast                ratio     SE  df asymp.LCL asymp.UCL null z.ratio p.value
# CFU0 / (CFU1e+05)       1.000 0.1602 Inf     0.682     1.468    1   0.003  1.0000
# CFU0 / (CFU1e+06)       0.156 0.0221 Inf     0.111     0.219    1 -13.126  <.0001
# (CFU1e+05) / (CFU1e+06) 0.156 0.0221 Inf     0.111     0.219    1 -13.097  <.0001
# 
# Confidence level used: 0.95 
# Conf-level adjustment: bonferroni method for 3 estimates 
# Intervals are back-transformed from the log scale 
# P value adjustment: bonferroni method for 3 tests 
# Tests are performed on the log scale 

car::Anova(me.model3)

# Analysis of Deviance Table (Type II tests)
# 
# Response: Surv(last.obs, censored)
# Df Chisq Pr(>Chisq)    
# CFU  2 252.3  < 2.2e-16 ***
# 	---
# 	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

######################
### Plots ###########
######################
# plot base model
ggsurv <- survminer::ggsurvplot(survfit(Surv(last.obs,censored)~CFU,data=exposure.aad),data=exposure.aad,
																conf.int = TRUE,censor = FALSE, palette="grey",
																legend.title = "Bacterial Load", legend.labs = c("Control", "Low", "High"))
ggsurv$plot + theme_classic()+
	theme(panel.grid.major.y = element_line(),
				panel.grid.minor.y = element_line(),
				axis.title = element_text(size = 14),
				axis.text = element_text(size = 14),
				legend.text = element_text(size = 14),
				legend.title = element_text(size = 14))+
	scale_y_continuous(limits = c(0,1))+
	scale_x_continuous(breaks = seq(0,50,5))+
	labs(x = "Exposure hours")

ggsave("grey_Cmac_adults_Vc_Pm_survplot.pdf", device = "pdf", path = ".",
			 dpi = 300, width = 150, height = 105, units = "mm")
