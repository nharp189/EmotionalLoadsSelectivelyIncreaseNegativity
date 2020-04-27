{library(lme4)
library(lmerTest)
library(tidyverse)
library(emmeans)
library(ggsignif)
library(papaja)
library(readxl)}

### set contrasts ###
options(contrasts = c("contr.sum","contr.poly"))


### Checking matrix characteristics ###
### import z score updated function ###
source("~/Documents/Nick-Grad/Neta_Lab/depletion_study/study2/Analyses/WorkingMemoryLoads/wilcox_test.R")
### code edits taken from : https://stats.stackexchange.com/questions/306841/z-score-on-wilcoxon-signed-ranks-test ###

stim.data <- readxl::read_xlsx("~/Documents/Nick-Grad/Neta_Lab/depletion_study/study2/Analyses/WorkingMemoryLoads/IAPS_Stim_List.xlsx")
shapiro.test(subset(stim.data, Condition == "POS")$Aro_Mn)
shapiro.test(subset(stim.data, Condition == "NEG")$Aro_Mn)
stim.t.test <- wilcox_test(subset(stim.data, Condition == "POS")$Aro_Mn, subset(stim.data, Condition == "NEG")$Aro_Mn)
stim.t.test$z_val
stim.t.test$p.value

### PRIMARY RESULTS ###
### Data from all trials (no memory probe exclusions) ###
data <- read.csv("~/Documents/Nick-Grad/Neta_Lab/depletion_study/study2/Analyses/WorkingMemoryLoads/Data/Cleaned_Data/Final_Data_NoExclusions.csv")

## drop first row
data <- data[, -1]

### check normality ###
shapiro.test(data$lo.emo.sur_rate)
shapiro.test(data$hi.emo.sur_rate)
shapiro.test(data$lo.neu.sur_rate)
shapiro.test(data$hi.neu.sur_rate)

shapiro.test(data$lo.emo.sur_n_MAD)
shapiro.test(data$lo.emo.sur_p_MAD)
shapiro.test(data$hi.emo.sur_n_MAD)
shapiro.test(data$hi.emo.sur_p_MAD)
shapiro.test(data$lo.neu.sur_n_MAD)
shapiro.test(data$lo.neu.sur_p_MAD)
shapiro.test(data$hi.neu.sur_n_MAD)
shapiro.test(data$hi.neu.sur_p_MAD)


### RATE ###
long <- gather(data, key = "Condition", value = "PerNeg",
               lo.emo.sur_rate, hi.emo.sur_rate,
               lo.neu.sur_rate, hi.neu.sur_rate)

level1 <- long[, c("subjID", "Condition", "PerNeg")]

### recode to numeric ###
level1$Dom <- ifelse(level1$Condition == "lo.emo.sur_rate", 1,
                     ifelse(level1$Condition == "hi.emo.sur_rate", 1, 0))
level1$Load <- ifelse(level1$Condition == "lo.emo.sur_rate", 0,
                      ifelse(level1$Condition == "lo.neu.sur_rate", 0, 1))

### make factors ###
level1$Dom <- dplyr::recode(level1$Dom,
                            "0" = "Non-emotional",
                            "1" = "Emotional")

level1$Load <- dplyr::recode(level1$Load,
                             "0" = "Low",
                             "1" = "High")
level1$subjID <- as.factor(level1$subjID)

### transform proportion to percent ###
level1$PerNeg.t <- level1$PerNeg * 100

rate.model <- (lmer(PerNeg.t ~ Dom * Load + (1 | subjID) + (1 | subjID:Dom) + (1 | subjID:Load), data = level1,
            REML = F))

### summary, ANOVA test, estimated marginal means, and 95% CIs
summary(rate.model)
anova(rate.model, type = "III")

### descriptive stats ###
mean(subset(level1, (level1$Dom == "Emotional"))$PerNeg.t, na.rm = T)
sd(subset(level1, (level1$Dom == "Emotional"))$PerNeg.t, na.rm = T)

mean(subset(level1, (level1$Dom == "Non-emotional"))$PerNeg.t, na.rm = T)
sd(subset(level1, (level1$Dom == "Non-emotional"))$PerNeg.t, na.rm = T)

mean(subset(level1, (level1$Load == "Low"))$PerNeg.t, na.rm = T)
sd(subset(level1, (level1$Load == "Low"))$PerNeg.t, na.rm = T)

mean(subset(level1, (level1$Load == "High"))$PerNeg.t, na.rm = T)
sd(subset(level1, (level1$Load == "High"))$PerNeg.t, na.rm = T)

### get domain estimates ###
emmeans(rate.model, pairwise ~ Dom, adjust = "none")
confint(emmeans(rate.model, pairwise ~ Dom, adjust = "none"))

### get load estimates ###
emmeans(rate.model, pairwise ~ Load, adjust = "none")
confint(emmeans(rate.model, pairwise ~ Load, adjust = "none"))


##############################################################
### maximum deviation analyses ###
long.md <- gather(data, key = "Condition", value = "MD",
                  lo.emo.sur_p_MAD, lo.emo.sur_n_MAD,
                  hi.emo.sur_p_MAD, hi.emo.sur_n_MAD,
                  lo.neu.sur_p_MAD, lo.neu.sur_n_MAD, 
                  hi.neu.sur_p_MAD, hi.neu.sur_n_MAD)

level1 <- long.md[, c("subjID", "Condition", "MD")]

level1$Dom <- ifelse(level1$Condition %in% c("lo.emo.sur_p_MAD",
                                             "hi.emo.sur_p_MAD",
                                             "lo.emo.sur_n_MAD",
                                             "hi.emo.sur_n_MAD"), "emotional", "non-emotional")

level1$Load <- ifelse(level1$Condition %in% c("lo.neu.sur_p_MAD",
                                              "lo.neu.sur_n_MAD",
                                              "lo.emo.sur_n_MAD",
                                              "lo.emo.sur_p_MAD"), "low", "high")

level1$Rate <- ifelse(level1$Condition %in% c("lo.emo.sur_p_MAD",
                                              "hi.emo.sur_p_MAD",
                                              "lo.neu.sur_p_MAD",
                                              "hi.neu.sur_p_MAD"), "positive", "negative")
level1$Rate <- as.factor(level1$Rate)
level1$Load <- as.factor(level1$Load)
level1$Dom <- as.factor(level1$Dom)
level1$subjID <- as.factor(level1$subjID)

md.mod <- lmer(MD ~  Rate * Load * Dom + (1 | subjID) + (1 | subjID:Rate) +
                  (1 | subjID:Load) + (1 | subjID:Dom), data = level1,
                REML = F)


### summary, ANOVA test, estimated marginal means, and 95% CIs
summary(md.mod)
anova(md.mod, type = "III")

### descriptive stats ###
## Main effect of Rate ###
mean(subset(level1, (level1$Rate == "positive"))$MD, na.rm = T)
sd(subset(level1, (level1$Rate == "positive"))$MD, na.rm = T)

mean(subset(level1, (level1$Rate == "negative"))$MD, na.rm = T)
sd(subset(level1, (level1$Rate == "negative"))$MD, na.rm = T)

### Rate x Load interaction ####
mean(subset(level1, (level1$Rate == "positive" & level1$Load == "low"))$MD, na.rm = T)
sd(subset(level1, (level1$Rate == "positive" & level1$Load == "low"))$MD, na.rm = T)

mean(subset(level1, (level1$Rate == "negative" & level1$Load == "low"))$MD, na.rm = T)
sd(subset(level1, (level1$Rate == "negative" & level1$Load == "low"))$MD, na.rm = T)

mean(subset(level1, (level1$Rate == "positive" & level1$Load == "high"))$MD, na.rm = T)
sd(subset(level1, (level1$Rate == "positive" & level1$Load == "high"))$MD, na.rm = T)

mean(subset(level1, (level1$Rate == "negative" & level1$Load == "high"))$MD, na.rm = T)
sd(subset(level1, (level1$Rate == "negative" & level1$Load == "high"))$MD, na.rm = T)

### main effect Domain ###
mean(subset(level1, (level1$Dom == "emotional"))$MD, na.rm = T)
sd(subset(level1, (level1$Dom == "emotional"))$MD, na.rm = T)

mean(subset(level1, (level1$Dom == "non-emotional"))$MD, na.rm = T)
sd(subset(level1, (level1$Dom == "non-emotional"))$MD, na.rm = T)

### Rate x Domain interaction ###
mean(subset(level1, (level1$Rate == "positive" & level1$Dom == "emotional"))$MD, na.rm = T)
sd(subset(level1, (level1$Rate == "positive" & level1$Dom == "emotional"))$MD, na.rm = T)

mean(subset(level1, (level1$Rate == "negative" & level1$Dom == "emotional"))$MD, na.rm = T)
sd(subset(level1, (level1$Rate == "negative" & level1$Dom == "emotional"))$MD, na.rm = T)

mean(subset(level1, (level1$Rate == "positive" & level1$Dom == "non-emotional"))$MD, na.rm = T)
sd(subset(level1, (level1$Rate == "positive" & level1$Dom == "non-emotional"))$MD, na.rm = T)

mean(subset(level1, (level1$Rate == "negative" & level1$Dom == "non-emotional"))$MD, na.rm = T)
sd(subset(level1, (level1$Rate == "negative" & level1$Dom == "non-emotional"))$MD, na.rm = T)


### main effect of rate ###
emmeans(md.mod, pairwise ~ Rate, adjust = "none")
confint(emmeans(md.mod, pairwise ~ Rate, adjust = "none"))

### load x rating ###
emmeans(md.mod, pairwise ~ Load:Rate, adjust = "none")
confint(emmeans(md.mod, pairwise ~ Load:Rate, adjust = "none"))

### main effect of domain ###
emmeans(md.mod, pairwise ~ Dom, adjust = "none")
confint(emmeans(md.mod, pairwise ~ Dom, adjust = "none"))

### Domain x rating ###
emmeans(md.mod, pairwise ~ Dom:Rate, adjust = "none")
confint(emmeans(md.mod, pairwise ~ Dom:Rate, adjust = "none"))

########################################


### test for differences on NEG vs POS probe accuracy 
iaps.norm <- readxl::read_excel("~/Documents/Nick-Grad/Neta_Lab/Stimuli/IAPS_ratings.xls")
stim.data <- readxl::read_xlsx("~/Documents/Nick-Grad/Neta_Lab/depletion_study/study2/Analyses/WorkingMemoryLoads/IAPS_Stim_List.xlsx")
### find the valence cutoffs for pos and neg in the 144 stim... 
min(stim.data$Val_Mn[which(stim.data$Condition == "POS")])
max(stim.data$Val_Mn[which(stim.data$Condition == "POS")])
min(stim.data$Val_Mn[which(stim.data$Condition == "NEG")])
max(stim.data$Val_Mn[which(stim.data$Condition == "NEG")])

### ok, now label the "no" images based off the valence cutoffs above... 
iaps.norm$CONDITION <- ifelse(iaps.norm$valmn >= 5.07, "POS",
                              ifelse(iaps.norm$valmn <= 4.32, "NEG",
                                     "NEU"))

### READ IN MT.data$DATA
MT.data$data$stim.recode2 <- str_sub(MT.data$data$stim, 22, 25)

MT.data$data$probeVal <- ifelse(MT.data$data$stim.recode2 %in% iaps.norm$IAPS,
                                iaps.norm$CONDITION, "")

names(MT.data$data)[names(MT.data$data)=="stim.recode2"] <- "IAPS"
MT.data$data <- merge(MT.data$data, iaps.norm[, c("IAPS", "CONDITION")], by = "IAPS")
MT.data$data$probeVal <- MT.data$data$CONDITION

MT.data$data <- subset(MT.data$data, MT.data$data$trialtype == "memory")
### fix 9635.2 discrepancy 
MT.data$data$probeVal <- ifelse(MT.data$data$IAPS == 9635, "NEG",
                                MT.data$data$probeVal)
MT.data$data$mem.cor <- as.numeric(MT.data$data$mem.cor)
MT.data.rating.table2 <- (ddply(MT.data$data, "subjID", summarise, 
                                emo.acc.Pos = mean(mem.cor[which(type == "EMO" & CONDITION == "POS")], na.rm = T),
                                emo.acc.Neg = mean(mem.cor[which(type == "EMO" & CONDITION == "NEG")], na.rm = T)))

t.test(MT.data.rating.table2$emo.acc.Pos,
       MT.data.rating.table2$emo.acc.Neg, paired = T)

### descriptives ###
mean(MT.data.rating.table2$emo.acc.Pos)
sd(MT.data.rating.table2$emo.acc.Pos)

mean(MT.data.rating.table2$emo.acc.Neg)
sd(MT.data.rating.table2$emo.acc.Neg)
