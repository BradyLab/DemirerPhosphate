#Figures for RSA of M82 and Penn for Paper and Posters
#Cutting previous work down to just bargraphs with jitter

library(lme4)
library(lmerTest)
library(lsmeans)
library(ggplot2)

setwd("C:/Users/Don/Desktop/Grad School/Brady Lab/Zinc Pi Paper/Root Hair Analysis/")

#This dataset is from the July 2016 phenotyping
finalroots <- read.csv("M82_Penn_Aug2018_RootHair_Processed.csv", header = TRUE)

#fit model for primary root length
lm1 <- lm(HairLenAve ~ Genotype*Phosphate, data = finalroots)
summary(lm1)
anova(lm1)
lsm.tmp <- lsmeans(lm1, ~ Genotype*Phosphate)
contrast(lsm.tmp, "pairwise", by = "Phosphate")
contrast(lsm.tmp, "trt.vs.ctrl", by = "Genotype")
contrast(lsm.tmp, "pairwise", interaction = T)

PL.results <- data.frame(
  Genotype=rep(levels(finalroots$Genotype),2),Phosphate=rep(levels(finalroots$Phosphate),1,each=nlevels(finalroots$Genotype)))

PL.results$Pred_PL_cm <- predict(lm1,PL.results,re.form=NA)
PL.results$sem.low <-PL.results$Pred_PL_cm - summary(lm1)$coefficients[,"Std. Error"]
PL.results$sem.high <- PL.results$Pred_PL_cm + summary(lm1)$coefficients[,"Std. Error"]
PL.results$t_value <- summary(lm1)$coefficients[,3]
PL.results$p_txt <- c("b","d","a","c")#summarized from the rsults of difflsmeans(lm1)
PL.results$y_stars <- PL.results$sem.high * 1.1
head(PL.results)

#set the reference
PL.results$Genotype <- relevel(PL.results$Genotype,ref="M82")
PL.results$Phosphate <- relevel(PL.results$Phosphate,ref="Sufficient")
#plot
pl <- ggplot(data=PL.results,aes(x=Phosphate,y=Pred_PL_cm))
pl <- pl + geom_bar(aes(color=Phosphate),stat="identity",position="dodge", fill = "#FFFFFF")
#pl <- pl + geom_jitter(data = finalroots,aes(x=Phosphate,y=PrimaryLen,color=Phosphate),size=0.5)#the original data
pl <- pl + geom_errorbar(aes(ymin = sem.low, ymax = sem.high,color=Phosphate), width = 0.5,position = "dodge") + scale_color_manual(values=c("#000000", "#CC33CC"))
pl <- pl + theme_classic()
pl <- pl + facet_grid(.~Genotype) + theme(text = element_text(size=16),axis.text = element_text(size = 16),axis.text.y = element_text(size = 16)) + ylab("Length (cm)") + geom_text(aes(y = y_stars,label=p_txt), position = position_dodge(width = 0.8))
pl <- pl + labs(title="Mature Root Hair Length") 
pl

#Lateral Root NUmber

lm2 <- lmer(NumberLat ~ Genotype*Phosphate + (1|Plate) , data = finalroots)
summary(lm2)
difflsmeans(lm2)

LR.results <- data.frame(
  Genotype=rep(levels(finalroots$Genotype),2),Phosphate=rep(levels(finalroots$Phosphate),1,each=nlevels(finalroots$Genotype)))

LR.results$Pred_LR_cm <- predict(lm2,LR.results,re.form=NA)
LR.results$sem.low <-LR.results$Pred_LR_cm - summary(lm2)$coefficients[,"Std. Error"]
LR.results$sem.high <- LR.results$Pred_LR_cm + summary(lm2)$coefficients[,"Std. Error"]
LR.results$t_value <- summary(lm1)$coefficients[,3]
LR.results$p_txt <- c("b","ab","a","c")#summarized from the rsults of difflsmeans(lm1)
LR.results$y_stars <- LR.results$sem.high * 1.1
head(LR.results)

LR.results$Genotype <- relevel(LR.results$Genotype,ref="M82")
LR.results$Phosphate <- relevel(LR.results$Phosphate,ref="Sufficient")
#plot
lr <- ggplot(data=LR.results,aes(x=Phosphate,y=Pred_LR_cm))
lr <- lr + geom_bar(aes(color=Phosphate),stat="identity",position="dodge", fill = "#FFFFFF")
#lr <- lr + geom_jitter(data = finalroots,aes(x=Phosphate,y=NumberLat,color=Phosphate),size=0.5)#the original data
lr <- lr + geom_errorbar(aes(ymin = sem.low, ymax = sem.high,color=Phosphate), width = 0.5,position = "dodge") + scale_color_manual(values=c("#000000", "#CC33CC"))
lr <- lr + theme_classic()
lr <- lr + facet_grid(.~Genotype) + theme(text = element_text(size=16),axis.text = element_text(size = 16),axis.text.y = element_text(size = 16)) + ylab("Number of Lateral Roots") + geom_text(aes(y = y_stars,label=p_txt), position = position_dodge(width = 0.8))
lr <- lr + labs(title="Lateral Root Number") 
lr

#Lateral Root Density calculated as Number lateral roots / primary root length

lm3 <- lmer(DensityLat ~ Genotype*Phosphate + (1|Plate) , data = finalroots)
summary(lm3)
difflsmeans(lm3)

LD.results <- data.frame(
  Genotype=rep(levels(finalroots$Genotype),2),Phosphate=rep(levels(finalroots$Phosphate),1,each=nlevels(finalroots$Genotype)))

LD.results$Pred_LD_cm <- predict(lm3,LD.results,re.form=NA)
LD.results$sem.low <-LD.results$Pred_LD_cm - summary(lm3)$coefficients[,"Std. Error"]
LD.results$sem.high <- LD.results$Pred_LD_cm + summary(lm3)$coefficients[,"Std. Error"]
LD.results$t_value <- summary(lm3)$coefficients[,3]
LD.results$p_txt <- c("a","c","a","b")#summarized from the rsults of difflsmeans(lm3)
LD.results$y_stars <- LD.results$sem.high * 1.1
head(LD.results)

LD.results$Genotype <- relevel(LD.results$Genotype,ref="M82")
LD.results$Phosphate <- relevel(LD.results$Phosphate,ref="Sufficient")
#plot
ld <- ggplot(data=LD.results,aes(x=Phosphate,y=Pred_LD_cm))
ld <- ld + geom_bar(aes(color=Phosphate),stat="identity",position="dodge", fill = "#FFFFFF")
#lr <- lr + geom_jitter(data = finalroots,aes(x=Phosphate,y=NumberLat,color=Phosphate),size=0.5)#the original data
ld <- ld + geom_errorbar(aes(ymin = sem.low, ymax = sem.high,color=Phosphate), width = 0.5,position = "dodge") + scale_color_manual(values=c("#000000", "#CC33CC"))
ld <- ld + theme_classic()
ld <- ld + facet_grid(.~Genotype) + theme(text = element_text(size=16),axis.text = element_text(size = 16),axis.text.y = element_text(size = 16)) + ylab("#Lat Root / 1' Root cm") + geom_text(aes(y = y_stars,label=p_txt), position = position_dodge(width = 0.8))
ld <- ld + labs(title="Lateral Root Density")
ld

# Total Root Length

lm4 <- lmer(TotalRootLen ~ Genotype*Phosphate + (1|Plate) , data = finalroots)
summary(lm4)
difflsmeans(lm4)


TRL.results <- data.frame(
  Genotype=rep(levels(finalroots$Genotype),2),Phosphate=rep(levels(finalroots$Phosphate),1,each=nlevels(finalroots$Genotype)))

TRL.results$Pred_TRL_cm <- predict(lm4,TRL.results,re.form=NA)
TRL.results$sem.low <-TRL.results$Pred_TRL_cm - summary(lm4)$coefficients[,"Std. Error"]
TRL.results$sem.high <- TRL.results$Pred_TRL_cm + summary(lm4)$coefficients[,"Std. Error"]
TRL.results$t_value <- summary(lm4)$coefficients[,3]
TRL.results$p_txt <- c("a","b","a","b")#summarized from the rsults of difflsmeans(lm4)
TRL.results$y_stars <- TRL.results$sem.high * 1.1
head(TRL.results)

TRL.results$Genotype <- relevel(TRL.results$Genotype,ref="M82")
TRL.results$Phosphate <- relevel(TRL.results$Phosphate,ref="Sufficient")
#plot
trl <- ggplot(data=TRL.results,aes(x=Phosphate,y=Pred_TRL_cm))
trl <- trl + geom_bar(aes(color=Phosphate),stat="identity",position="dodge", fill = "#FFFFFF")
#lr <- lr + geom_jitter(data = finalroots,aes(x=Phosphate,y=NumberLat,color=Phosphate),size=0.5)#the original data
trl <- trl + geom_errorbar(aes(ymin = sem.low, ymax = sem.high,color=Phosphate), width = 0.5,position = "dodge") + scale_color_manual(values=c("#000000", "#CC33CC"))
trl <- trl + theme_classic()
trl <- trl + facet_grid(.~Genotype) + theme(text = element_text(size=16),axis.text = element_text(size = 16),axis.text.y = element_text(size = 16)) + ylab("Length (cm)") + geom_text(aes(y = y_stars,label=p_txt), position = position_dodge(width = 0.8))
trl <- trl + labs(title="Total Root Length")
trl

# Total Lateral Length

lm5 <- lmer(TotalLatLen ~ Genotype*Phosphate + (1|Plate) , data = finalroots)
summary(lm5)
difflsmeans(lm5)


LRL.results <- data.frame(
  Genotype=rep(levels(finalroots$Genotype),2),Phosphate=rep(levels(finalroots$Phosphate),1,each=nlevels(finalroots$Genotype)))

LRL.results$Pred_LRL_cm <- predict(lm5,LRL.results,re.form=NA)
LRL.results$sem.low <-LRL.results$Pred_LRL_cm - summary(lm5)$coefficients[,"Std. Error"]
LRL.results$sem.high <- LRL.results$Pred_LRL_cm + summary(lm5)$coefficients[,"Std. Error"]
LRL.results$t_value <- summary(lm5)$coefficients[,3]
LRL.results$p_txt <- c("b","a","a","a")#summarized from the rsults of difflsmeans(lm5)
LRL.results$y_stars <- LRL.results$sem.high * 1.1
head(LRL.results)

LRL.results$Genotype <- relevel(LRL.results$Genotype,ref="M82")
LRL.results$Phosphate <- relevel(LRL.results$Phosphate,ref="Sufficient")
#plot
LRL <- ggplot(data=LRL.results,aes(x=Phosphate,y=Pred_LRL_cm))
LRL <- LRL + geom_bar(aes(color=Phosphate),stat="identity",position="dodge", fill = "#FFFFFF")
#lr <- lr + geom_jitter(data = finalroots,aes(x=Phosphate,y=NumberLat,color=Phosphate),size=0.5)#the original data
LRL <- LRL + geom_errorbar(aes(ymin = sem.low, ymax = sem.high,color=Phosphate), width = 0.5,position = "dodge") + scale_color_manual(values=c("#000000", "#CC33CC"))
LRL <- LRL + theme_classic()
LRL <- LRL + facet_grid(.~Genotype) + theme(text = element_text(size=16),axis.text = element_text(size = 16),axis.text.y = element_text(size = 16)) + ylab("Length (cm)") + geom_text(aes(y = y_stars,label=p_txt), position = position_dodge(width = 0.8))
LRL <- LRL + labs(title="Total Lateral Root Length")
LRL
