library(ggplot2)
library(lme4)
library(moments)
library(MASS)
library(bbmle) 
library(broom.mixed)
library(caret)
library(multcomp)

###########################################
#### import data and data prep ############
###########################################

data.in <- read.csv("confocal_JIPipe_Cmac_Hoechst_algae_size.csv",stringsAsFactors = TRUE)
names(data.in)
data.in$Algae <- as.numeric(data.in$Algae)
data.in$Hoechst <- as.numeric(data.in$Hoechst)
data.in$Exposure.days <- as.numeric(data.in$Exposure.days)
data.in$Replicate_fixation_day <- as.factor(data.in$Replicate_fixation_day)
data.in$Condition <- as.factor(data.in$Condition)
data.in$Individual <- as.factor(data.in$Individual)
data.in$orientation <-  as.factor(data.in$orientation)
data.in$Length.mm <- as.numeric(data.in$Length.mm)

# get ratio algal cells / animal cells
data.in$Ratio.algae.hoechst <- data.in$Algae / data.in$Hoechst


#############################################################
#### stats - linear model - continuous variable #############
#############################################################
### only 2 days data ###


# only do model for 2 days samples
data.2 <- data.in[data.in$Exposure.days == 2,]

data.2 <- droplevels(data.2)

lm.2d.1 <-lm(Ratio.algae.hoechst ~ Condition * Length.mm * orientation , data = data.2)
res.2d <- resid(lm.2d.1)
shapiro.test(res.2d)

# Shapiro-Wilk normality test
# 
# data:  res.2d
# W = 0.94726, p-value = 0.02232

# no transformation needed, residuals are already normal

skewness(res.2d)
# 0.8978783
# skewed to the left but no substantial departure from normality (<2.1 see Kim 2019)

kurtosis(res.2d)
# 4.172668
# leptokurtic but no substantial departure from normality (<7 see Kim 2019)

hist(res,breaks=20)
# skewness clearly seen, but looks okay

caret::findLinearCombos(model.matrix(lm.2d.1))
# $remove
# [1] 9 15

colnames(model.matrix(lm.2d.1)) [c(9, 15)]


# try removing interactions
lm.2d.2 <- lm(Ratio.algae.hoechst ~ Condition + Length.mm + orientation , data = data.2)

anova(lm.2d.1,lm.2d.2)
# p = 0.6433 > 0.01 =? non-significant, keep simpler model 

car::Anova(lm.2d.2)

# try removing orientation
lm.2d.3 <- lm(Ratio.algae.hoechst ~ Condition + Length.mm  , data = data.2)

anova(lm.2d.2,lm.2d.3)
# p = 0.4869 > 0.01 -> non-significant, keep simpler model

# try removing length
lm.2d.4 <- lm(Ratio.algae.hoechst ~ Condition  , data = data.2)

anova(lm.2d.3,lm.2d.4)
# p = 0.1405 > 0.01 -> non significant, keep simpler model

# compare with null model

lm.2d.null <- lm(Ratio.algae.hoechst ~ 1  , data = data.2)

anova(lm.2d.4, lm.2d.null)
# p =  0.006961 < 0.01 -> significantly different, keep more complex model with condition

summary(lm.2d.4)
# Call:
#   lm(formula = Ratio.algae.hoechst ~ Condition, data = data.2)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.1941 -0.9271 -0.5835  0.8292  5.8180 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      3.1678     0.3817   8.299  6.7e-11 ***
#   Condition1e+05  -0.7388     0.5482  -1.348  0.18394    
# Condition1e+06  -1.7328     0.5254  -3.298  0.00182 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.574 on 49 degrees of freedom
# (23 observations deleted due to missingness)
# Multiple R-squared:  0.1835,	Adjusted R-squared:  0.1502 
# F-statistic: 5.507 on 2 and 49 DF,  p-value: 0.006961

car::Anova(lm.2d.4)

# Anova Table (Type II tests)
#
# Response: Ratio.algae.hoechst
# Sum Sq Df F value   Pr(>F)   
# Condition  27.282  2  5.5069 0.006961 **
#   Residuals 121.376 49                    


ggplot(data = data.2,
       aes(y = Ratio.algae.hoechst, x = Condition))+
  geom_violin(fill = "grey80", colour = "grey80")+
  geom_point(size = 2)+
  theme_minimal()+
  theme(panel.grid.major.y = element_line(),
      panel.grid.minor.y = element_line(),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14))+
  labs(x = "\nBacterial Load", y = "Ratio algae/animal cells\n")+
  scale_x_discrete(labels= c("Control","Low dose","High dose"))


ggsave("ratio_algae_hoechst_2days.pdf", device = "pdf", path = ".",
       dpi = 300, width = 140 , height = 105 , units = "mm")

### post hoc comparisons 

# Tukey
summary(glht(lm.2d.4, linfct = mcp(Condition = "Tukey")))

# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: Tukey Contrasts
# 
# 
# Fit: lm(formula = Ratio.algae.hoechst ~ Condition, data = data.2)
# 
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)   
# 1e+05 - 0 == 0      -0.7388     0.5482  -1.348  0.37598   
# 1e+06 - 0 == 0      -1.7328     0.5254  -3.298  0.00514 **
#   1e+06 - 1e+05 == 0  -0.9940     0.5340  -1.861  0.16080   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)

# # T test with Bonferroni correction
# pairwise.t.test(data.2$Ratio.algae.hoechst, data.2$Condition, p.adjust.method = "bonferroni")
# 
# # Pairwise comparisons using t tests with pooled SD 
# # 
# # data:  data.2$Ratio.algae.hoechst and data.2$Condition 
# # 
# # 0      1e+05 
# # 1e+05 0.5518 -     
# #   1e+06 0.0055 0.2061
# # 
# # P value adjustment method: bonferroni 

#############################
### check animal lengths ####
############################

ggplot(data = data.2, aes(x=Condition,y=Length.mm))+
  geom_violin(fill = "grey80", colour = "grey80")+
  geom_point(size = 2)+
  theme_minimal()+
  theme(panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))+
  labs(x = "\nBacterial Load", y = "Animal length\n")+
  scale_x_discrete(labels= c("Control","Low dose","High dose"))

ggsave("Cmac_length_2days.pdf", device = "pdf", path = ".",
       dpi = 300, width = 140 , height = 105 , units = "mm")


lm.length <-lm(Length.mm ~ Condition , data = data.2)
res.l <- resid(lm.length)
shapiro.test(res.l)

# Shapiro-Wilk normality test
# 
# data:  res.l
# W = 0.98377, p-value = 0.6956
# not non-normal

# compare with null model
null.l.m <- lm(Length.mm ~ 1 , data = data.2)

anova(lm.length,null.l.m)
# p = 0.001231 < 0.05

# equivalent to
car::Anova(lm.lenght)

# Anova Table (Type II tests)
# 
# Response: Length.mm
# Sum Sq Df F value   Pr(>F)   
# Condition 17.349  2  7.7059 0.001231 **
#   Residuals 55.160 49                    

# Tukey post hoc
summary(glht(lm.length, linfct = mcp(Condition = "Tukey")))

# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: Tukey Contrasts
# 
# 
# Fit: lm(formula = Length.mm ~ Condition, data = data.2)
# 
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)    
# 1e+05 - 0 == 0      -0.9691     0.3696  -2.622   0.0306 *  
#   1e+06 - 0 == 0      -1.3625     0.3542  -3.847   <0.001 ***
#   1e+06 - 1e+05 == 0  -0.3934     0.3600  -1.093   0.5229    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)


