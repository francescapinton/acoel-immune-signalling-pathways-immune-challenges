library(ggplot2)
library(survminer)
library(ggstance)
library(car)
library(MASS)

# read in the exposure data
data.init <- read.csv2(file = "Cmac-longterm-input-data.csv",stringsAsFactors = T)

names(data.init)
levels(as.factor(data.init$Date.start))
levels(as.factor(data.init$days.after.feeding))
levels(as.factor(data.init$days.exposure))
levels(as.factor(data.init$CFU))
levels(as.factor(data.init$Dead))
levels(as.factor(data.init$Offspring))
length(levels(as.factor(data.init$Individual)))

# remove Notes column (empty)
exposure.data <- data.init[1:(length(data.init)-1)]
# make sure data are ordered 
exposure.data <- exposure.data[order(exposure.data$Individual, exposure.data$days.exposure), ]

###############################
#### Cox survival analyses ####
###############################
library(coxme)
library(survival)
library(ggplot2)
library(survminer)
library(emmeans)

#data prep
# one row per individual, recording the first date at which it was seen dead (min)
# since Dead value is included, there will be two rows per individuals, one with Dead 1 and one with Dead 0 (unless they were always alive or dead from the beginning)
last.obs.min <- aggregate(exposure.data$days.exposure,by=list(id = exposure.data$Individual, Bacteria = exposure.data$Bacteria, CFU = exposure.data$CFU, Date.start = exposure.data$Date.start,Dead=exposure.data$Dead),min,na.rm=TRUE)
# remove all the alive status rows
last.obs.dead <- last.obs.min[last.obs.min$Dead == 1,]
names(last.obs.dead)[names(last.obs.dead) == "x"] <- "last.obs"
length(levels(last.obs.dead$id))
last.obs.dead <- droplevels(last.obs.dead)
length(levels(last.obs.dead$id))
# length is for checking - it is normal that some levels are dropped
# (the ids of individuals that were alive at the end of the experiment)

# same but with max to get the last timepoint where they were seen alive (or dead)
last.obs.max <- aggregate(exposure.data$days.exposure,by=list(id = exposure.data$Individual, Bacteria = exposure.data$Bacteria, CFU = exposure.data$CFU, Date.start = exposure.data$Date.start,Dead=exposure.data$Dead),max)
last.obs.alive <- last.obs.max[last.obs.max$Dead == 0,]
names(last.obs.alive)[names(last.obs.alive) == "x"] <- "last.obs"
length(levels(last.obs.alive$id))
last.obs.alive <- droplevels(last.obs.alive)
length(levels(last.obs.alive$id))
# here only individuals that were dead at 0 d should be removed
# but they were already removed manually because "absent-from-start" in the notes column

# combine dead and alive dataframes AND only keep only the first instance (last observation from the dead df) if both present 
# (both are present only if the individual died, otherwise three will only be one occurrence in the alive df)
last.obs.both <- rbind(last.obs.dead,last.obs.alive)
last.obs <- last.obs.both[match(unique(last.obs.both$id),last.obs.both$id),]
last.obs <- droplevels(last.obs)
length(levels(last.obs$id))

# CFUs as factors
last.obs$CFU <- as.factor(last.obs$CFU)


### Cox mixed effects model
# survival object
inf.surv <- Surv(last.obs$last.obs,last.obs$Dead)

# null model
null.model <- survival::coxph(inf.surv~1)
summary(null.model)

# base model
# time of last observation with censored status (survival object) as a function of CFU
# CFU is considered a discrete variable (factor), not continuous
base.model <- survival::coxph(inf.surv~last.obs$CFU)
summary(base.model)

# test models
anova(null.model,base.model)
# significantly different p = 8.947e-06 -> reject simpler model

# mixed effect model
me.model <- coxme::coxme(inf.surv~ CFU + (1 | Date.start /id),data = last.obs)
summary(me.model)

# simplify model: remove id
me.model1 <- coxme::coxme(inf.surv~ CFU + (1 | Date.start ),data = last.obs)
summary(me.model1)

anova(me.model,me.model1)
# p=0.9841 -> not significantly different -> keep simpler model

# simplify further: remove random effect term --> base model
anova(me.model1,base.model)
# p=1.592e-08 -> significantly different
# keep model with random effect: significantly better than base model

# minimal adequate model is me.model1
summary(me.model1)

# Mixed effects coxme model
# Formula: inf.surv ~ CFU + (1 | Date.start) 
# Data: last.obs 
# 
# events, n = 166, 216
# 
# Random effects:
# 	group  variable        sd  variance
# 1 Date.start Intercept 0.6633183 0.4399912
# Chisq  df         p   AIC   BIC
# Integrated loglik 55.19 3.0 6.266e-12 49.19 39.85
# Penalized loglik 63.17 3.9 5.387e-13 55.36 43.22
# 
# Fixed effects:
# 	coef exp(coef) se(coef)    z        p
# CFU100000  0.2473    1.2806   0.1989 1.24    0.214
# CFU1000000 0.9975    2.7115   0.1975 5.05 4.41e-07

car::Anova(me.model1) 
# Analysis of Deviance Table (Type II tests)
# 
# Response: inf.surv
# Df  Chisq Pr(>Chisq)    
# CFU  2 28.646   6.02e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## post hoc comparisons

me.model1.posthoc <-  emmeans(me.model1, specs = pairwise ~ CFU, type = "response",adjust = "bonferroni")
summary(me.model1.posthoc,infer = TRUE)

# $emmeans
# CFU     hazard     SE  df asymp.LCL asymp.UCL null z.ratio p.value
# 0        0.660 0.0766 Inf     0.526     0.829    1  -3.577  0.0003
# 100000   0.846 0.0944 Inf     0.679     1.053    1  -1.501  0.1333
# 1000000  1.791 0.1980 Inf     1.441     2.225    1   5.256  <.0001
# 
# Confidence level used: 0.95 
# Intervals are back-transformed from the log scale 
# Tests are performed on the log scale 
# 
# $contrasts
# contrast               ratio     SE  df asymp.LCL asymp.UCL null z.ratio p.value
# CFU0 / CFU100000       0.781 0.1550 Inf     0.485     1.257    1  -1.244  0.6409
# CFU0 / CFU1000000      0.369 0.0728 Inf     0.230     0.592    1  -5.050  <.0001
# CFU100000 / CFU1000000 0.472 0.0897 Inf     0.300     0.744    1  -3.951  0.0002
# 
# Confidence level used: 0.95 
# Conf-level adjustment: bonferroni method for 3 estimates 
# Intervals are back-transformed from the log scale 
# P value adjustment: bonferroni method for 3 tests 
# Tests are performed on the log scale 

# plot model without random effects
ggsurv <- survminer::ggsurvplot(survfit(inf.surv~CFU,data=last.obs),data=last.obs,
																conf.int = TRUE,censor = FALSE, palette="grey")
ggsurv$plot + theme_classic()+
	theme(panel.grid.major.y = element_line(),
				panel.grid.minor.y = element_line(),
				axis.title = element_text(size = 14),
				axis.text = element_text(size = 14))+
	scale_x_continuous(breaks = seq(0,85,5))+
	labs(x = "Exposure days")+
	theme(legend.position = "none")

ggsave("grey_Cmac_adults_Vc_survival_survplot_longterm.pdf", device = "pdf", path = ".",
			 dpi = 300, width = 150, height = 105, units = "mm")


