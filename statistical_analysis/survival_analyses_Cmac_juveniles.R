library(ggplot2)
library(survminer)
library(ggstance)

# read in the exposure data
prel.exposure.data <- read.csv(file = "20240411_Vcoralliilyticus_juveniles.csv",stringsAsFactors = T)
prel.exposure.data <- rbind(prel.exposure.data,read.csv(file = "20240414_Vcoralliilyticus_juveniles.csv",stringsAsFactors = T))
prel.exposure.data <- rbind(prel.exposure.data,read.csv(file = "20240419_Vcoralliilyticus_juveniles.csv",stringsAsFactors = T))
prel.exposure.data <- rbind(prel.exposure.data,read.csv(file = "20240424_Vcoralliilyticus_juveniles.csv",stringsAsFactors = T))

names(prel.exposure.data)
levels(as.factor(prel.exposure.data$Date.start))
# substitute NA in the Notes column with empty string so that these rows are not remove while subsetting
prel.exposure.data[is.na(prel.exposure.data$Notes),]$Notes <- ""

length(levels(prel.exposure.data$Individual))
# remove juveniles not present at 0d of exposure
exposure.data <- prel.exposure.data[prel.exposure.data$Notes != "absent-from-start", ]
exposure.data <- droplevels(exposure.data)
# check if it dropped more individuals than expected
length(levels(exposure.data$Individual))

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

######################################
### all days observed (up to 9 or 10)
#######################################

### data prep: for each individual: age at death or last seen and status (censored or not)

# sum dead values for each individual
# inf.status <- aggregate(exposure.data$Dead,by = list(id = exposure.data$Individual, Bacteria = exposure.data$Bacteria, CFU = exposure.data$CFU, Date.start = exposure.data$Date.start),sum)
# # rename the sum column censored and if it is greater than 1 (i.e. the individual is dead at some point), make it 1 (censored status = dead)
# names(inf.status)[names(inf.status) == "x"] <- "censored"
# inf.status[inf.status$censored > 0, ]$censored <-  1

# create a dataframe with only status dead or alive at the end of recording period
# not needed anymore
# inf.last.obs <- aggregate(exposure.data$days.exposure,by = list(id = exposure.data$Individual, Bacteria = exposure.data$Bacteria, CFU = exposure.data$CFU, Date.start = exposure.data$Date.start),last.obs.f)

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

# consider bacterial dose as discrete (factor)
last.obs$CFU <- as.factor(last.obs$CFU)


### Cox mixed effects model
inf.surv <- Surv(last.obs$last.obs,last.obs$Dead)

# null model
null.model <- survival::coxph(inf.surv~1)
summary(null.model)

# basic Cox's proportional hazard model:
# time of last observation with censored status (survival object) as a function of CFU
# CFU is considered a discrete variable (factor), not continuous
base.model <- survival::coxph(inf.surv~last.obs$CFU)
summary(base.model)

#### build and test models
anova(null.model,base.model)
# p = 0.06428
# base.model NOT significantly better than null model

# plot base model
ggsurv <- survminer::ggsurvplot(survfit(inf.surv~CFU,data=last.obs),data=last.obs,
																conf.int = TRUE,censor = FALSE, palette="grey",
																legend.title = "Bacterial Load", legend.labs = c("Control", "Low", "High"))
ggsurv$plot + theme_classic()+
	theme(panel.grid.major.y = element_line(),
				panel.grid.minor.y = element_line(),
				axis.title = element_text(size = 14),
				axis.text = element_text(size = 14),
				legend.text = element_text(size = 14),
				legend.title = element_text(size = 14))+
	scale_x_continuous(breaks = seq(0,10,1))+
	labs(x = "Exposure days")+
	theme(legend.position = "none")
	
ggsave("grey_Cmac_juveniles_Vc_survival_survplot_10d_all_CFUs.pdf", device = "pdf", path = ".",
			 dpi = 300, width = 150, height = 105, units = "mm")

# plot null model
ggsurv <- survminer::ggsurvplot(survfit(inf.surv~1,data=last.obs),data=last.obs,
																conf.int = TRUE,censor = FALSE, color="black")
ggsurv$plot + theme_classic()+
	theme(panel.grid.major.y = element_line(),
				panel.grid.minor.y = element_line(),
				axis.title = element_text(size = 14),
				axis.text = element_text(size = 14))+
	scale_x_continuous(breaks = seq(0,10,1))+
	labs(x = "Exposure days")

ggsave("grey_Cmac_juveniles_Vc_survival_survplot_10d_one_only.pdf", device = "pdf", path = ".",
			 dpi = 300, width = 150, height = 105, units = "mm")

#############################################################
# analyses for the first 2 days only

# remove all days of exposure after 2
exposure.data.2d <- exposure.data[exposure.data$days.exposure < 3,]
exposure.data.2d <- droplevels(exposure.data.2d)
levels(as.factor(exposure.data.2d$days.exposure))
length(levels(exposure.data.2d$Individual))


### data prep: for each individual: age at death or last seen and status (censored or not)

# one row per individual, recording the first date at which it was seen dead (min)
# since Dead value is included, there will be two rows per individuals, one wiht Dead 1 and one with Dead 0 (unless they were always alive or dead from the beginning)
last.obs.min.2d <- aggregate(exposure.data.2d$days.exposure,by=list(id = exposure.data.2d$Individual, Bacteria = exposure.data.2d$Bacteria, CFU = exposure.data.2d$CFU, Date.start = exposure.data.2d$Date.start,Dead=exposure.data.2d$Dead),min,na.rm=TRUE)
# remove all the alive status rows
last.obs.dead.2d <- last.obs.min.2d[last.obs.min.2d$Dead == 1,]
names(last.obs.dead.2d)[names(last.obs.dead.2d) == "x"] <- "last.obs"
length(levels(last.obs.dead.2d$id))
last.obs.dead.2d <- droplevels(last.obs.dead.2d)
length(levels(last.obs.dead.2d$id))
# length is for checking - it is normal that some levels are dropped
# (the ids of individuals that were alive at the end of the experiment)

# same but with max to get the last timepoint where they were seen alive (or dead)
last.obs.max.2d <- aggregate(exposure.data.2d$days.exposure,by=list(id = exposure.data.2d$Individual, Bacteria = exposure.data.2d$Bacteria, CFU = exposure.data.2d$CFU, Date.start = exposure.data.2d$Date.start,Dead=exposure.data.2d$Dead),max)
last.obs.alive.2d <- last.obs.max.2d[last.obs.max.2d$Dead == 0,]
names(last.obs.alive.2d)[names(last.obs.alive.2d) == "x"] <- "last.obs"
length(levels(last.obs.alive.2d$id))
last.obs.alive.2d <- droplevels(last.obs.alive.2d)
length(levels(last.obs.alive.2d$id))
# here only individuals that were dead at 0 d should be removed
# but they were already removed manually because "absent-from-start" in the notes column

# combine dead and alive dataframes AND only keep only the first instance (last observation from the dead df) if both present 
# (both are present only if the individual died, otherwise three will only be one occurrence in the alive df)
last.obs.both.2d <- rbind(last.obs.dead.2d,last.obs.alive.2d)
last.obs.2d <- last.obs.both.2d[match(unique(last.obs.both.2d$id),last.obs.both.2d$id),]
last.obs.2d <- droplevels(last.obs.2d)
length(levels(last.obs.2d$id))

# bacterial dose discrete (factor)
last.obs.2d$CFU <- as.factor(last.obs.2d$CFU)

## create survival object
inf.surv.2d <- Surv(last.obs.2d$last.obs,last.obs.2d$Dead)


### Cox mixed effects model
# null model
null.model.2d <- survival::coxph(inf.surv.2d~1)
summary(null.model.2d)

# basic Cox's proportional hazard model:
# time of last observation with censored status (survival object) as a function of CFU
# CFU is considered a discrete variable (factor), not continuous
base.model.2d <- survival::coxph(inf.surv.2d~last.obs.2d$CFU)
summary(base.model.2d)

#### build and test models
anova(null.model.2d,base.model.2d)
# base model not significantly better than null model 
# p = 0.07153

# mixed effects model
# Fixed effects: CFU
# Random effects: experiment/replicate (Date.start), Individual
me.model.2d <- coxme::coxme(inf.surv.2d~CFU + (1 | Date.start / id ),data = last.obs.2d)
summary(me.model.2d)

anova(base.model.2d,me.model.2d)
# mixed effect model not significantly better than base model p= 0.1745

# test assumptions
# 1) common baseline hazard by design since all animals come from the same aquarium and same culture conditions
# 2) proportional hazard assumption with alpha = 0.001
# i.e. residuals are independent of time (time invariant covariates and regression coefficients) - correct??
(me.model.2d.test <- survival::cox.zph(me.model.2d))
# p > 0.001 -> proportional hazards
# if alpha = 0.05 -> non-proportional hazards
# https://stats.stackexchange.com/questions/499735/violated-non-proportional-hazards-cox-regression-model-of-time-dependent-covar
# https://stats.stackexchange.com/questions/359970/checking-the-proportional-hazard-assumption/400981#400981
# interpret them as weighted average of hazard ratios - they are still different and change over time
# not trusting the value itself
survminer::ggcoxzph(me.model.2d.test,df = 2)


# test for normality
res.2d <- resid(coxph(inf.surv.2d ~ last.obs.2d$CFU),data=last.obs.2d)
hist(res.2d,breaks=50)
shapiro.test(res.2d)
# not normal residual distribution - Not necessary for cox model

# simplify the me model
# remove id
me.model.2d.2 <- coxme::coxme(inf.surv.2d~CFU + (1 | Date.start ),data = last.obs.2d)
summary(me.model.2d.2)

anova(me.model.2d,me.model.2d.2)
# non significantly different -> accept the simpler model me.model.2d.2

# remove the fixed coefficient CFU
me.model.2d.3 <- coxme::coxme(inf.surv.2d~1 + (1 | Date.start ),data = last.obs.2d)
summary(me.model.2d.3)

anova(me.model.2d.2,me.model.2d.3)
# not significantly different --> accept the simpler model me.model.2d.3

# compare me.model.2d.3 with base model (more simplified version)
anova(base.model.2d,me.model.2d.3)
# p = 0.1985 -> no significant difference
# keep base model
# but base model is not significantly different than null model
anova(base.model.2d,null.model.2d)

# keep null model
summary(null.model.2d)
# Call:  survival::coxph(formula = inf.surv.2d ~ 1)


##################
#### Plots #######
#################

# plot model with CFU as explanatory variable
ggsurv.2d <- survminer::ggsurvplot(survfit(inf.surv.2d~CFU,data=last.obs.2d),data=last.obs.2d,
																conf.int = TRUE,censor = FALSE, palette="grey",
																legend.title = "Bacterial Load", legend.labs = c("Control", "Low", "High"))
ggsurv.2d$plot + theme_classic()+
	theme(panel.grid.major.y = element_line(),
				panel.grid.minor.y = element_line(),
				axis.title = element_text(size = 14),
				axis.text = element_text(size = 14),
				legend.text = element_text(size = 14),
				legend.title = element_text(size = 14))+
	scale_x_continuous(breaks = seq(0,2,1))+
	scale_y_continuous(limits = c(0.50,1))+
	labs(x = "Exposure days")+
	theme(legend.position = "none")

ggsave("grey_Cmac_juveniles_Vc_survival_survplot_2d_all_CFUs.pdf", device = "pdf", path = ".",
			 dpi = 300, width = 150, height = 105, units = "mm")


# plot null model
ggsurv.2d <- survminer::ggsurvplot(survfit(inf.surv.2d~1,data=last.obs.2d),data=last.obs.2d,
																	 conf.int = TRUE,censor = FALSE, color = "black")

ggsurv.2d$plot + theme_classic()+
	theme(panel.grid.major.y = element_line(),
				panel.grid.minor.y = element_line(),
				axis.title = element_text(size = 14),
				axis.text = element_text(size = 14))+
	scale_x_continuous(breaks = seq(0,2,1))+
	scale_y_continuous(limits = c(0.50,1))+
	labs(x = "Exposure days")

ggsave("grey_Cmac_juveniles_Vc_survival_survplot_2d_one_only.pdf", device = "pdf", path = ".",
			 dpi = 300, width = 150, height = 105, units = "mm")
