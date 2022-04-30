#Daily temperature variationmagnifies the toxicity of amixture consisting of a chemical pesticide and a biopesticide in a vector mosquito
#Vienna Delnat, Tam Tran, Lizanne Janssens and Robby Stoks
#Science of the Total Environment (2019)
#R code tested on 16/03/2022

#####Packages#####

install.packages("afex")
install.packages("lme4")
install.packages("car")     
install.packages("lsmeans")
install.packages("effects")     

library(afex)
library(lme4)
library(car)     
library(lsmeans)
library(effects)     

sessionInfo()

##Set working directory to source file
#RStudio -> Session -> Set Working Directory...-> To Source File Location

#####Datasets#####

#Dataset voor TotalMortality
dataBIN=read.csv("Delnat-et-al_DTV-larvae-mortality-binomial.csv", sep=",", na.strings=c(""))
##Set correct data types 
Factors <- c("DTV", "Pesticide", "Bti", "CPF", "day", "number", "replicate", "id")
dataBIN[Factors] <- do.call(cbind.data.frame, lapply(dataBIN[Factors], as.factor))
##Set levels in factor
dataBIN$Pesticide=factor(dataBIN$Pesticide,levels=c("Control","Bti","Chlorpyrifos","Mixture")) 
#Subset
dataBIN=subset(dataBIN, Mortality2d!="NA")
str(dataBIN) 

#Dataset voor WingLength and DevelopmentTime
#CombinedDTVsolved --> with new wing lengths (no animals used with AdultDate=31,32,33)
dataADULT=read.csv("./Delnat-et-al_DTV-adults-sublethal.csv", sep=",", na.strings=c(""))
##Set correct data types 
Factors <- c("DTV", "Pesticide", "Bti", "CPF", "day", "number", "id", "Jar", "Epje", "StartDate", "AdultDate", "Dead", "Sex", "SideWing")
dataADULT[Factors] <- do.call(cbind.data.frame, lapply(dataADULT[Factors], as.factor))
##Set levels in factor
dataADULT$Pesticide=factor(dataADULT$Pesticide,levels=c("CTR","Bti65","Bti70","CPF","Mix65","Mix70")) 
str(dataADULT) 


#####Total mortality - Larvae#####

#Interaction in model --> use set_sum_contrasts() and type=3 in Anova
set_sum_contrasts()

#Generalized linear mixed models with a binomial error structure and the logit link
#Correct for pseudoreplication (replicate = Vial)
glmerTotalMortality=glmer(TotalMortality ~ DTV*CPF*Bti + (1|replicate), data=dataBIN, na.action=na.omit, family=binomial(link=logit),
                          control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e5)))
Anova(glmerTotalMortality, type="III") 

#Posthoc test - contrast analysis with fdr correction
pairs(lsmeans(glmerTotalMortality, ~CPF*Bti|DTV, adjust="fdr")) 
pairs(lsmeans(glmerTotalMortality, ~DTV|CPF*Bti, adjust="fdr")) 

#Quick effects plot - not used in manuscript
plot(effect(mod=glmerTotalMortality, term="DTV*CPF*Bti"))
#Emmeans and standard errors for figure
TotalMortalityPlotData <- summary(lsmeans(glmerTotalMortality, ~ DTV*CPF*Bti, type = "response"))

#Assumption - Dispersion parameter
glmTotalMortality=glm(TotalMortality ~ DTV*CPF*Bti, data=dataBIN, na.action=na.omit, family=quasibinomial(link=logit))
summary(glmTotalMortality) 


#####Development time - adult females#####

#Subset
dataADULTfemale=subset(dataADULT, Sex=="female")

#Interaction in model --> use set_sum_contrasts() and type=3 in Anova
set_sum_contrasts()

#Correct for pseudoreplication (Jar = Vial)
lmerDevelopmentTime=lmer(DevelopmentTimeEgg ~ DTV*CPF*Bti + (1|Jar), data=dataADULTfemale, na.action=na.omit)
Anova(lmerDevelopmentTime, type="III", white.adjust=TRUE) 

#Posthoc test - contrast analysis with fdr correction
pairs(lsmeans(lmerDevelopmentTime, ~CPF*Bti|DTV, adjust="fdr")) 
pairs(lsmeans(lmerDevelopmentTime, ~DTV|CPF*Bti, adjust="fdr")) 

#Quick effects plot - not used in manuscript
plot(effect(mod=lmerDevelopmentTime, term="DTV*CPF*Bti"))
#Emmeans and standard errors for figure
DevelopmentTimeEggPlotData <- summary(emmeans(lmerDevelopmentTime, ~ DTV*CPF*Bti, type = "response"))

#Assumption - Normality of residuals
shapiro.test(resid(lmerDevelopmentTime))                  
hist(resid(lmerDevelopmentTime))      
#Assumption - Homogeneity of variance
leveneTest(DevelopmentTime ~ DTV*CPF*Bti, data = dataADULTfemale)
#Thumb of rule - if minimum and maximum variance do not differ more than a factor 5 - assumption still met
aggregate(DevelopmentTime ~ DTV*CPF*Bti, data = dataADULTfemale, var) 

#Outliers and influential observations
outlierTest(lmerDevelopmentTime)
cd=cooks.distance(lmerDevelopmentTime); which(cd>1)
influenceIndexPlot(lmerDevelopmentTime, vars = c("studentized", "Bonf"))


#####Wing length - adult females#####

#Interaction in model --> use set_sum_contrasts() and type=3 in Anova
set_sum_contrasts()

#General linear mixed model with a normal error structure and the identity link
#Correct for pseudoreplication (Jar = Vial)
lmerWingLength=lmer(WingLength ~ DTV*CPF*Bti + (1|Jar), data=dataADULT, na.action=na.omit)
Anova(lmerWingLength, type="III") 

#Posthoc test - contrast analysis with fdr correction
pairs(lsmeans(lmerWingLength, ~CPF*Bti|DTV, adjust="fdr")) 
pairs(lsmeans(lmerWingLength, ~DTV|CPF*Bti, adjust="fdr")) 

#Quick effects plot - not used in manuscript
plot(effect(mod=lmerWingLength, term="DTV*CPF*Bti"))
#Emmeans and standard errors for figure
WingLengthPlotData <- summary(emmeans(lmerWingLength, ~ DTV*CPF*Bti, type = "response"))

#Assumption - Normality of residuals
shapiro.test(resid(lmerWingLength))                  
hist(resid(lmerWingLength))      
#Assumption - Homogeneity of variance
leveneTest(WingLength ~ DTV*CPF*Bti, data = dataADULTfemale)
#Thumb of rule - if minimum and maximum variance do not differ more than a factor 5 - assumption still met
aggregate(WingLength ~ DTV*CPF*Bti, data = dataADULTfemale, var) 

#Outliers and influential observations
outlierTest(lmerWingLength)
cd=cooks.distance(lmerWingLength); which(cd>1)
influenceIndexPlot(lmerWingLength, vars = c("studentized", "Bonf"))


######Save Rdata######
save.image(file="Chapter3_Rdata_20220316_NotPublished.Rdata")
