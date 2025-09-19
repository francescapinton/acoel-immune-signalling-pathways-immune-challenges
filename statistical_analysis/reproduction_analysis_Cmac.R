library(ggplot2)
library(moments)
library(lme4)
library(magrittr)
library(car)
library(dplyr)
library(tidyr)
library(DHARMa)
library(glmmADMB)
library(glmmTMB)
library(bbmle) 
library(emmeans)

#### data prep ####

# read in the exposure data
exposure.data <- read.csv(file = "Cmac-2d-input-data.csv",stringsAsFactors = T)
names(exposure.data)

# consider only surely NON-CONTAMINATED replicates and NOT partial replicate without 10^6 CFU
exposure.data <- exposure.data[exposure.data$Date.start != "2023-01-25",]
exposure.data <- exposure.data[exposure.data$Date.start != "2023-05-24",]
exposure.data <- droplevels(exposure.data)

exposure.data1 <- exposure.data

#NA in cell suspension or bud column considered as 0
exposure.data1["bud"][is.na(exposure.data1["bud"])] <- 0
exposure.data1$bud <- as.factor(exposure.data1$bud)
exposure.data1["cells.suspension"][is.na(exposure.data1["cells.suspension"])] <- 0
exposure.data1$cells.suspension <- as.factor(exposure.data1$cells.suspension)

# categorical variables
exposure.data1$bud <- as.factor(exposure.data1$bud)

# take max num of (corrected) produced offspring for each individual
# corrected offspring = offspring - offspring at 0h (i.e. produced before start of exposure)
offsprings <- aggregate(exposure.data1$corr.Offspring,
                        by = list(individual = exposure.data1$Individual, CFU = exposure.data1$CFU, bud = exposure.data1$bud, Date.start = exposure.data1$Date.start),
                        FUN = max)
names(offsprings)[names(offsprings) == "x"] <- "num.os"


os.data.2d <- exposure.data1[,c("Date.start","CFU","Individual","Dead","corr.Offspring","hours.pi","bud")]
os.data.2d$CFU <- as.factor(os.data.2d$CFU)

# remove rows with NA for corr.offspring (i.e. individuals dead at that timepoint)
os.data.2d <-  os.data.2d[!is.na(os.data.2d$corr.Offspring),]

# keep only data for 24h and 48h
os.data.2d <- os.data.2d[os.data.2d$hours.pi %in% c(24,48),]
levels(as.factor(os.data.2d$hours.pi))

# days instead of hours
os.data.2d$days.exposure <- 1
os.data.2d[os.data.2d$hours.pi == 48,]$days.exposure <- 2

# calculate difference for days.exposure = 2
# order data correctly
os.data.2d <- os.data.2d[order(os.data.2d$Individual, os.data.2d$days.exposure),]
os.data.2d$os.produced <- ave(os.data.2d$corr.Offspring, os.data.2d$Individual, FUN = function(x) c(x[1], diff(x)))

levels(as.factor(os.data.2d$os.produced))

# remove all animals dead with offspring (should have not been recorded)
os.data.2d <- os.data.2d[os.data.2d$Dead == 0,]

# offspring production as continuous response variable
os.data.2d$os.produced <- as.numeric(os.data.2d$os.produced)

# discrete explanatory variable (too few levels)
os.data.2d$days.exposure <- as.factor(os.data.2d$days.exposure)
os.data.2d$CFU <- as.factor(os.data.2d$CFU)
os.data.2d$bud <- as.factor(os.data.2d$bud)

# not considered in the model, but also discrete
os.data.2d$Dead <- as.factor(os.data.2d$Dead)

os.data.2d <- droplevels(os.data.2d)

summary(os.data.2d)

# split measurements for day 1 and day 2 because repeated measurements create temporal pseudoreplication
# Crawley 2013 p 683-684
os.data.d1 <- os.data.2d[os.data.2d$days.exposure == 1,]
os.data.d2 <- os.data.2d[os.data.2d$days.exposure == 2,]

### test assumptions ###

# fit lm
os.lm.d1 <- lm(os.produced ~ CFU * bud, data = os.data.d1)
os.lm.d2<- lm(os.produced ~ CFU * bud, data = os.data.d2)

# check for normality
res.d1<-resid(os.lm.d1)
hist(res.d1,breaks=30)
shapiro.test(res.d1)

res.d2<-resid(os.lm.d2)
hist(res.d2,breaks=30)
shapiro.test(res.d2)

# Shapiro-Wilk normality test
# for both p-value < 2.2e-16
# rejected H0 hypothesis of normal distribution

skewness(res.d1) # measure of asymmetry see Kim 2013
# 2.027813 < 2.1
kurtosis(res.d1) # measure of peakedness see Kim 2013
# 7.460499 < 7.5
skewness(res.d2) # measure of asymmetry see Kim 2013
# 2.745184 > 2.1
kurtosis(res.d2) # measure of peakedness see Kim 2013
# 10.68096 > 7.5
# d1 barely inside the thresholds,
# d2 substantially different from normality
# given Shapiro-Wilk consider both of them non-normal

# plot of residuals and fitted values
plot_lm.d1 = os.lm.d1
plot_lm.d1 %>% plot(which = 1)
plot_lm.d2 = os.lm.d2
plot_lm.d2 %>% plot(which = 1)

# try Box-Cox transformation to obtain a normal distribution
# try Box-Cox transformation shifting response variable to make it positive (adding +1)
library(MASS)
bc.d1 <- boxcox(lm(os.produced +1 ~ CFU * bud, data = os.data.d1),lambda = seq(-10, 10, 0.1))
bc.d2 <- boxcox(lm(os.produced +1 ~ CFU * bud, data = os.data.d2),lambda = seq(-10, 10, 0.1))

# find optimal lambda
bc.d1$x[which.max(bc$y)] # -6.9 -> super low, not meaningful
bc.d2$x[which.max(bc$y)] # -6.9 -> super low, not meaningful

### check for overdispersion and zero inflation ####
# index of dispersion
sum(residuals(os.lm.d1, type="pearson")^2) / df.residual(os.lm.d1) # 0.1842716
sum(residuals(os.lm.d2, type="pearson")^2) / df.residual(os.lm.d2) # 0.134462

# under-dispersed for both
# see https://en.wikipedia.org/wiki/Index_of_dispersion

# # check for zero inflation
# sim_res <- simulateResiduals(os.lm)
# testZeroInflation(sim_res)

testZeroInflation(os.lm.d1)
testZeroInflation(os.lm.d2)
# ratioObsSim = Inf, p-value < 2.2e-16
# they are zero inflated - but should be done after to check, not before to determine error distribution
# see https://search.r-project.org/CRAN/refmans/DHARMa/html/testZeroInflation.html

# it makes sense biologically that it is zero inflated - many animals did not reproduce (some are not mature and some have not reproduced in this time window)
# needs a zero inflated model?

# choose procedure
# non-normal data -> no lm, no lmer
# needs random effects -> no lm, no glm
# not overdispersed but zero inflated-> use glmmADMB or glmmTMB, rather than glmer 


#### determine family distribution ####
### day 1
# poisson
m.d1.p <- glmmTMB(os.produced ~ CFU * bud + (1 | Date.start), data = os.data.d1,
									ziformula=~0,family=poisson())
summary(m.d1.p)

# poisson zero-inflated
m.d1.p.zi <- glmmTMB(os.produced ~ CFU * bud + (1 | Date.start), data = os.data.d1,
									ziformula=~1,family=poisson())
summary(m.d1.p.zi)

# negative binomial 1
m.d1.nb1 <- glmmTMB(os.produced ~ CFU * bud + (1 | Date.start), data = os.data.d1,
									ziformula=~0,family=nbinom1())
summary(m.d1.nb1)

# negative binomial 1 zero-inflated
m.d1.nb1.zi <- glmmTMB(os.produced ~ CFU * bud + (1 | Date.start), data = os.data.d1,
										 ziformula=~1,family=nbinom1())
#   Model convergence problem; non-positive-definite Hessian matrix.
#  Model convergence problem; false convergence (8)
summary(m.d1.nb1.zi)


# negative binomial 2
m.d1.nb2 <- glmmTMB(os.produced ~ CFU * bud + (1 | Date.start), data = os.data.d1,
										ziformula=~0,family=nbinom2())
summary(m.d1.nb2)

# negative binomial 2 zero-inflated
m.d1.nb2.zi <- glmmTMB(os.produced ~ CFU * bud + (1 | Date.start), data = os.data.d1,
											 ziformula=~1,family=nbinom2())
summary(m.d1.nb2.zi)

AICtab(m.d1.p,m.d1.p.zi,m.d1.nb1,m.d1.nb1.zi,m.d1.nb2,m.d1.nb2.zi)

#            dAIC df
# m.d1.p       0   7 
# m.d1.p.zi    2   8 
# m.d1.nb2     2   8 
# m.d1.nb1     2   8 
# m.d1.nb2.zi  4   9 
# m.d1.nb1.zi NA   9 

# poisson is the best, but small delta 
# -> substantial support for models with delta < 2 and strong support < 4 (Burnham and Anderson 2004)

# try glmer with poisson just to check
# control --> maximum likelihood for comparison
# glmmTMB uses ML as default (REML = FALSE by default)
glmer.d1.p <- glmer(os.produced ~ CFU * bud + (1 | Date.start), data = os.data.d1,
									family=poisson(), control = glmerControl(optimizer = "bobyqa"))
AIC(m.d1.p,glmer.d1.p)
# exactly the same AIC

### day 2
# poisson
m.d2.p <- glmmTMB(os.produced ~ CFU * bud + (1 | Date.start), data = os.data.d2,
									ziformula=~0,family=poisson())
summary(m.d2.p)

# poisson zero-inflated
m.d2.p.zi <- glmmTMB(os.produced ~ CFU * bud + (1 | Date.start), data = os.data.d2,
										 ziformula=~1,family=poisson())
summary(m.d2.p.zi)

# negative binomial 1
m.d2.nb1 <- glmmTMB(os.produced ~ CFU * bud + (1 | Date.start), data = os.data.d2,
										ziformula=~0,family=nbinom1())
summary(m.d2.nb1)

# negative binomial 1 zero-inflated
m.d2.nb1.zi <- glmmTMB(os.produced ~ CFU * bud + (1 | Date.start), data = os.data.d2,
											 ziformula=~1,family=nbinom1())
summary(m.d2.nb1.zi)


# negative binomial 2
m.d2.nb2 <- glmmTMB(os.produced ~ CFU * bud + (1 | Date.start), data = os.data.d2,
										ziformula=~0,family=nbinom2())
summary(m.d2.nb2)

# negative binomial 2 zero-inflated
m.d2.nb2.zi <- glmmTMB(os.produced ~ CFU * bud + (1 | Date.start), data = os.data.d2,
											 ziformula=~1,family=nbinom2())
#  Model convergence problem; non-positive-definite Hessian matrix.
summary(m.d2.nb2.zi)

AICtab(m.d2.p,m.d2.p.zi,m.d2.nb1,m.d2.nb1.zi,m.d2.nb2,m.d2.nb2.zi)

#            dAIC df
# m.d2.p      0.0  7 
# m.d2.nb1    1.5  8 
# m.d2.p.zi   2.0  8 
# m.d2.nb2    2.0  8 
# m.d2.nb1.zi 3.5  9 
# m.d2.nb2.zi  NA  9 

# very similar to d1
# poisson distribution is the best

#### simplify model ####
# comparisons between models with different effects are okay because ML estimation is used by default
### day 1
summary(m.d1.p)
car::Anova(m.d1.p)

# remove interaction (non-significant)
m.d1.p.1 <- glmmTMB(os.produced ~ CFU + bud + (1 | Date.start), data = os.data.d1,
									ziformula=~0,family=poisson())
car::Anova(m.d1.p.1)

anova(m.d1.p, m.d1.p.1)
# p = 0.8075 - non significantly different from each other
# -> keep the simpler model (m.d1.p.1)

# remove CFU (non significant term)
m.d1.p.2 <- glmmTMB(os.produced ~ bud + (1 | Date.start), data = os.data.d1,
										ziformula=~0,family=poisson())
car::Anova(m.d1.p.2)

anova(m.d1.p.1, m.d1.p.2)
# p = 0.2994 -> keep simpler model (m.d1.p.2)

# remove random effects
m.d1.p.3 <- glmmTMB(os.produced ~ bud , data = os.data.d1,
										ziformula=~0,family=poisson())
car::Anova(m.d1.p.3)

anova(m.d1.p.2, m.d1.p.3)
# p = 0.02016 < 0.05 - significant 
# keep more complex model with random effects

# best model m.d1.p.2
summary(m.d1.p.2)

# Family: poisson  ( log )
# Formula:          os.produced ~ bud + (1 | Date.start)
# Data: os.data.d1
# 
# AIC      BIC   logLik deviance df.resid 
# 858.2    872.3   -426.1    852.2      810 
# 
# Random effects:
# 	
# 	Conditional model:
# 	Groups     Name        Variance Std.Dev.
# Date.start (Intercept) 0.06351  0.252   
# Number of obs: 813, groups:  Date.start, 4
# 
# Conditional model:
# 	Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.7658     0.1557 -11.339  < 2e-16 ***
# 	bud1          0.9026     0.1817   4.968 6.75e-07 ***
# 	---
# 	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

car::Anova(m.d1.p.2)

# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: os.produced
# Chisq Df Pr(>Chisq)    
# bud 24.684  1  6.752e-07 ***
# 	---
# 	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# day 2
summary(m.d2.p)
car::Anova(m.d2.p)

# remove non-significant interaction CFU * bud
m.d2.p.1 <-  glmmTMB(os.produced ~ CFU + bud + (1 | Date.start), data = os.data.d2,
										 ziformula=~0,family=poisson())
summary(m.d2.p.1)
anova(m.d2.p,m.d2.p.1)
# p = 0.9712 -> not significantly different - keep simpler model m.d2.p.1
car::Anova(m.d2.p.1)

# remove CFU (non-significant term)
m.d2.p.2 <-  glmmTMB(os.produced ~ bud + (1 | Date.start), data = os.data.d2,
										 ziformula=~0,family=poisson())
summary(m.d2.p.2)

anova(m.d2.p.1,m.d2.p.2)
# p = 0.1317 -> not significantly different - keep simpler model m.d2.p.2

#remove random effects
m.d2.p.3 <-  glmmTMB(os.produced ~ bud, data = os.data.d2,
										 ziformula=~0,family=poisson())
anova(m.d2.p.3,m.d2.p.2)
# 0.1163 -> not significantly different
# keep simpler model without random effects
summary(m.d2.p.3)

# Family: poisson  ( log )
# Formula:          os.produced ~ bud
# Data: os.data.d2
# 
# AIC      BIC   logLik deviance df.resid 
# 532.6    541.6   -264.3    528.6      667 
# 
# 
# Conditional model:
# 	Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -2.2212     0.1250 -17.770  < 2e-16 ***
# 	bud1          1.0298     0.2394   4.303 1.69e-05 ***
# 	---
# 	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

car::Anova(m.d2.p.3)

# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: os.produced
# Chisq Df Pr(>Chisq)    
# bud 18.512  1  1.688e-05 ***
# 	---
# 	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

