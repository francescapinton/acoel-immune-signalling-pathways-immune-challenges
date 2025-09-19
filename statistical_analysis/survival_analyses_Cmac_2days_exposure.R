# read in the exposure data
exposure.data <- read.csv(file = "Cmac-2d-input-data.csv",stringsAsFactors = T)
names(exposure.data)

#NA in cell suspension or bud column considered as 0
exposure.data["bud"][is.na(exposure.data["bud"])] <- 0
exposure.data$bud <- as.factor(exposure.data$bud)
exposure.data["cells.suspension"][is.na(exposure.data["cells.suspension"])] <- 0
exposure.data$cells.suspension <- as.factor(exposure.data$cells.suspension)

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
exposure.aad[d.sum$x == 4,]$last.obs <- 4
exposure.aad[d.sum$x == 5,]$last.obs <- 2
exposure.aad[d.sum$x == 6,]$last.obs <- 1

hist(exposure.aad$last.obs)

### Cox mixed effects model
# null model
null.model <- survival::coxph(Surv(last.obs,censored)~1,data=exposure.aad)
summary(null.model)

# basic Cox's proportional hazard model:
# time of last observation with censored status (survival object) as a function of CFU
# CFU is considered a discrete variable (factor), not continuous
base.model <- survival::coxph(Surv(last.obs,censored)~CFU,data=exposure.aad)
summary(base.model)

#### build and test models
anova(null.model,base.model)
# base.model significantly better than null model

# Fixed effects: CFU, cell.suspension
# Random effects: experiment/replicate (Date.start), Individual

me.model <- coxme::coxme(Surv(last.obs,censored)~CFU*initial.damage + (1 | Date.start / id ),data = exposure.aad)
summary(me.model)

#### test assumptions
# 1) common baseline hazard by design since all animals come from the same aquarium and same culture conditions
# 2) proportional hazard assumption with alpha = 0.001
# i.e. residuals are independent of time (time invariant covariates and regression coefficients) - correct??
(me.model.test <- survival::cox.zph(me.model))
# p > 0.001 -> proportional hazards
# if alpha = 0.05 -> non-proportional hazards
# https://stats.stackexchange.com/questions/499735/violated-non-proportional-hazards-cox-regression-model-of-time-dependent-covar
# https://stats.stackexchange.com/questions/359970/checking-the-proportional-hazard-assumption/400981#400981
# interpret them as weighted average of hazard ratios - they are still different and change over time
# not trusting the value itself
survminer::ggcoxzph(me.model.test,df = 2)
plot(me.model.test)


# test for normality
res <- resid(coxph(Surv(exposure.aad$last.obs, exposure.aad$censored) ~ exposure.aad$CFU + exposure.aad$initial.damage),data=exposure.aad)
hist(res,breaks=50)
shapiro.test(res)
# not normal residual distribution - SHOULD RESIDUALS BE NORMAL FOR COX PROP HAZARD MODEL?
# Bewick 2013: "In Cox's model no assumption is made about the probability distribution of the hazard."


### simplify the model
# remove the (non-significant, p > 0.001 in all effects with interaction) interaction between CFU and initial.damage
me.model1 <- coxme::coxme(Surv(last.obs,censored)~CFU + initial.damage + (1 | Date.start / id ),data = exposure.aad)
summary(me.model1)

anova(me.model,me.model1)
# p> 0.1 -> accept the simpler model me.model1

# no id as nested random effect
me.model2 <- coxme::coxme(Surv(last.obs,censored)~CFU + initial.damage + (1 | Date.start),data = exposure.aad)
summary(me.model2)

anova(me.model1,me.model2)
# p > 0.6 -> accept simpler model without id

# removed cells in suspension factor
me.model3 <- coxme::coxme(Surv(last.obs,censored)~CFU + (1 | Date.start ),data = exposure.aad)
summary(me.model3)
anova(me.model2,me.model3)
# p > 0.05 -> accept simpler model without initial damage

# remove random effects completely
fe.model <- survival::coxph(Surv(last.obs,censored)~CFU,data = exposure.aad)
summary(fe.model)

anova(me.model3,fe.model)
# p << 0.001 -> random effect of experiment Start date is important

# compare me.model3 and base model
anova(base.model,me.model3)
# p << 0.001 -> me.model3 kept

#### results ####
# minimal adequate model is me.model3
summary(me.model3)

# Mixed effects coxme model
# Formula: Surv(last.obs, censored) ~ CFU + (1 | Date.start) 
# Data: exposure.aad 
# 
# events, n = 192, 863
# 
# Random effects:
# 	group  variable        sd  variance
# 1 Date.start Intercept 0.6924304 0.4794598
# Chisq   df p   AIC   BIC
# Integrated loglik 128.3 3.00 0 122.3 112.5
# Penalized loglik 140.2 4.85 0 130.5 114.7
# 
# Fixed effects:
# 	coef exp(coef) se(coef)    z        p
# CFU1e+05 0.8017    2.2293   0.2268 3.53 0.000409
# CFU1e+06 1.5361    4.6465   0.2104 7.30 2.85e-13

car::Anova(me.model3)

# Analysis of Deviance Table (Type II tests)
# 
# Response: Surv(last.obs, censored)
# Df  Chisq Pr(>Chisq)    
# CFU  2 59.742  1.065e-13 ***

###################################
### post hoc comparisons ##########
###################################

# interaction-style plot: visualize the impact of the different factors on the response estimated marginal means (means computed from the model)
emmip(ref_grid(me.model2),initial.damage~CFU,style = "factor",type = "response")
# not valid anymore for me.model2

me.model3.posthoc <-  emmeans(me.model3, specs = pairwise ~ CFU, type = "response",adjust = "bonferroni")
summary(me.model3.posthoc,infer = TRUE)
 
# $emmeans
# CFU   response     SE  df asymp.LCL asymp.UCL null z.ratio p.value
# 0        0.459 0.0621 Inf     0.352     0.599    1  -5.755  <.0001
# 1e+05    1.024 0.1140 Inf     0.823     1.274    1   0.209  0.8345
# 1e+06    2.133 0.2140 Inf     1.753     2.597    1   7.554  <.0001
# 
# Confidence level used: 0.95 
# Intervals are back-transformed from the log scale 
# Tests are performed on the log scale 
# 
# $contrasts
# contrast                ratio     SE  df asymp.LCL asymp.UCL null z.ratio p.value
# CFU0 / (CFU1e+05)       0.449 0.1020 Inf     0.261     0.772    1  -3.534  0.0012
# CFU0 / (CFU1e+06)       0.215 0.0453 Inf     0.130     0.356    1  -7.301  <.0001
# (CFU1e+05) / (CFU1e+06) 0.480 0.0783 Inf     0.325     0.709    1  -4.500  <.0001
# 
# Confidence level used: 0.95 
# Conf-level adjustment: bonferroni method for 3 estimates 
# Intervals are back-transformed from the log scale 
# P value adjustment: bonferroni method for 3 tests 
# Tests are performed on the log scale 


library(bbmle)
AICtab(me.model3,base.model)

############################
### plot survival curves ###
############################

# not possible to plot coxme model, so plot coxph only with fixed effects and  and without random effects
# require(survival)
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
	scale_y_continuous(limits = c(0.50,1))+
	scale_x_continuous(breaks = seq(0,50,5))+
	labs(x = "Exposure hours")+
	theme(legend.position = "none")

ggsave("grey_Cmac_adults_Vc_survival_survplot_legend.pdf", device = "pdf", path = ".",
			 dpi = 300, width = 150, height = 105, units = "mm")



# scale_color_manual(values = c("grey80","grey40","black"))


# Plot Kaplan-Meier survival curve
km_fit <- survfit(Surv(last.obs,censored) ~ CFU, data= exposure.aad)
plot(km_fit)

ggsurvplot(km_fit,risk.table = TRUE, 
					 pval = TRUE, conf.int = TRUE,
					 legend.title = "CFUs", legend.labs = c("0","10^5","10^6"))


##########################
### Hofstenia miamia ####
#########################


exposure.data.Hmia <- read.csv(file = "Hmia-input-data.csv",stringsAsFactors = T)
exposure.data.Hmia <- exposure.data.Hmia[,!(names(exposure.data.Hmia) %in% c("bud","cells.suspension"))]
names(exposure.data.Hmia)[names(exposure.data.Hmia) == "Offspring"] <- "Eggs"
names(exposure.data.Hmia)[names(exposure.data.Hmia) == "corr.Offspring"] <- "corr.Eggs"
names(exposure.data.Hmia)

# exclude very high dose and very low dose (only pilots)
exposure.data.Hmia <- exposure.data.Hmia[exposure.data.Hmia$CFU != 10000,]
exposure.data.Hmia <- exposure.data.Hmia[exposure.data.Hmia$CFU != 25000000,]

### data prep: for each individual: age at death or last seen and status (censored or not)
# sum dead values for each individual
d.sum.hmia <- aggregate(exposure.data.Hmia$Dead,by = list(id = exposure.data.Hmia$Individual, Bacteria = exposure.data.Hmia$Bacteria, CFU = exposure.data.Hmia$CFU, Date.start = exposure.data.Hmia$Date.start),sum)
# create a new dataframe with one row per individual
exposure.aad.hmia <- data.frame(d.sum.hmia[,c(1,2,3,4,5)])

# give censored status 0 to all the ones seen alive at 48h (sum of dead observations 0), 1 to all the others
exposure.aad.hmia$censored <-  1
exposure.aad.hmia[d.sum.hmia$x == 0,]$censored <- 0
# consider CFU as a categorical variable
exposure.aad.hmia$CFU <- as.factor(exposure.aad.hmia$CFU)

# add last observation (death or alive for the ones at 48h)
exposure.aad.hmia$last.obs <- NA
exposure.aad.hmia[d.sum.hmia$x == 0 | d.sum.hmia$x == 1,]$last.obs <- 48
exposure.aad.hmia[d.sum.hmia$x == 2,]$last.obs <- 24
exposure.aad.hmia[d.sum.hmia$x == 3,]$last.obs <- 6
exposure.aad.hmia[d.sum.hmia$x == 4,]$last.obs <- 4
exposure.aad.hmia[d.sum.hmia$x == 5,]$last.obs <- 2
exposure.aad.hmia[d.sum.hmia$x == 6,]$last.obs <- 1

### Cox mixed effects model
# null model
null.model.hmia <- survival::coxph(Surv(last.obs,censored)~1,data=exposure.aad.hmia)
summary(null.model.hmia)

ggsurv.Hmia <- survminer::ggsurvplot(survfit(Surv(last.obs,censored)~CFU,data=exposure.aad.hmia),data=exposure.aad.hmia,
																conf.int = TRUE,censor = FALSE, color ="#333333",
																legend.title = "Bacterial Load", legend.labs = c("Control", "Low", "High"))
ggsurv.Hmia$plot + theme_classic()+
	theme(panel.grid.major.y = element_line(),
				panel.grid.minor.y = element_line(),
				axis.title = element_text(size = 14),
				axis.text = element_text(size = 14),
				legend.text = element_text(size = 14),
				legend.title = element_text(size = 14))+
	scale_x_continuous(breaks = seq(0,50,5))+
	scale_y_continuous(limits = c(0.50,1))+
	labs(x = "Exposure hours")


ggsave("Hmia_Vc_survival_survplot_high.pdf", device = "pdf", path = ".",
			 dpi = 300, width = 140 , height = 105 , units = "mm")

ggsurv.Hmia <- survminer::ggsurvplot(survfit(Surv(last.obs,censored)~CFU,data=exposure.aad.hmia),data=exposure.aad.hmia,
																		 conf.int = TRUE,censor = FALSE, color ="#cccccc",
																		 legend.title = "Bacterial Load", legend.labs = c("Control", "Low", "High"))
ggsurv.Hmia$plot + theme_classic()+
	theme(panel.grid.major.y = element_line(),
				panel.grid.minor.y = element_line(),
				axis.title = element_text(size = 14),
				axis.text = element_text(size = 14),
				legend.text = element_text(size = 14),
				legend.title = element_text(size = 14))+
	scale_x_continuous(breaks = seq(0,50,5))+
	scale_y_continuous(limits = c(0.50,1))+
	labs(x = "Exposure hours")


ggsave("Hmia_Vc_survival_survplot_control.pdf", device = "pdf", path = ".",
			 dpi = 300, width = 140 , height = 105 , units = "mm")


ggsurv.Hmia <- survminer::ggsurvplot(survfit(Surv(last.obs,censored)~CFU,data=exposure.aad.hmia),data=exposure.aad.hmia,
																		 conf.int = TRUE,censor = FALSE, color ="#989898",
																		 legend.title = "Bacterial Load", legend.labs = c("Control", "Low", "High"))
ggsurv.Hmia$plot + theme_classic()+
	theme(panel.grid.major.y = element_line(),
				panel.grid.minor.y = element_line(),
				axis.title = element_text(size = 14),
				axis.text = element_text(size = 14),
				legend.text = element_text(size = 14),
				legend.title = element_text(size = 14))+
	scale_x_continuous(breaks = seq(0,50,5))+
	scale_y_continuous(limits = c(0.50,1))+
	labs(x = "Exposure hours")


ggsave("Hmia_Vc_survival_survplot_low.pdf", device = "pdf", path = ".",
			 dpi = 300, width = 140 , height = 105 , units = "mm")

