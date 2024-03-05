#CCL Target Trial 
#Does surgical management for cranial cruciate ligament (CCL) rupture cause a difference in lameness and analgesia outcomes compared to non-surgical management?
#IPW used to adjust for confounding and IPCW to adjust for censoring
#Outcomes as if all dogs treated and all dogs not treated predicted
#E-values used as sensitivity analysis
#Missing indicator approach for missing data
#Data including adult-only cases and <1 month between diagnosis and surgery 

#Packages to use
library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape)
library(finalfit)
library(knitr)
library(epiDisplay)
library(sjPlot)
library(pROC)
library(broom)
library(ipw)
library(forcats)
library(survey)
library(broom)
library(tableone)
library(Matching)
library(reshape2)
library(boot)
library(epitools)
library(foreign)

#Import raw data
CCL_dat <- read.csv("C:/PhD/230913_CCL_Target_Trial_Open_Access.csv",header=T,stringsAsFactors = T)


#Group Breeds
CCL_dat$VetCompassBreed <- fct_other(CCL_dat$VetCompassBreed_Cleaned, keep=c('Crossbreed', 'Labrador Retriever', 'Jack Russell Terrier', 'West Highland White Terrier',
                                                                             'Staffordshire Bull Terrier', 'Golden Retriever', 'English Springer Spaniel',
                                                                             'Cocker Spaniel', 'Yorkshire Terrier', 'Bichon Frise', 'Rottweiler'
                                                                              ),other_level="Other")

#Categorise bodyweight (-100 = bodyweight not recorded)
CCL_dat$Bodyweightcat <- cut(CCL_dat$Bodyweight..kg., breaks=c(-101, 0, 15, 70))

CCL_dat$Bodyweight_cat <- factor(CCL_dat$Bodyweightcat, levels=c("(0,15]","(15,70]","(-101,0]"))

#Convert necessary variables to factors
CCL_dat$Bodyweight..kg._missing <- factor(CCL_dat$Bodyweight..kg._missing)
CCL_dat$Overweight.obese_missing <- factor(CCL_dat$Overweight.obese_missing)

#estimate inverse probability weights (IPWs) (using package). Covariates to include derived from a DAG.
weight <- ipwpoint(
  exposure = Management.of.CCL.rupture_01,
  family = "binomial",
  link = "logit",
  numerator = ~ 1,
  denominator = ~ Bodyweightcat + Bodyweight..kg._missing + Age.at.diagnosis..years. 
  + Overweight.obese + Overweight.obese_missing + Orthopaedic.comorbidity.at.diagnosis.collapsed + Number.of.non.orthopaedic.comorbidities.at.diagnosis_grouped
  + Insurance.status + DataSiloName + Neuter_Status + VetCompassBreed,
  data = CCL_dat)

summary(weight$ipw.weights)

#Plot weights
ipwplot(weights = weight$ipw.weights, logscale = FALSE, 
        main = "Stabilized weights", xlim = c(0, 8))

#Outcome models

#Short-term lameness
#Generate inverse probability of censoring weights (IPCW) for lameness at 3 months
denom.cens <- glm(Lameness_3mo_cens ~ Bodyweightcat + Bodyweight..kg._missing + Age.at.diagnosis..years. 
                  + Overweight.obese + Overweight.obese_missing + Orthopaedic.comorbidity.at.diagnosis.collapsed + Number.of.non.orthopaedic.comorbidities.at.diagnosis_grouped
                  + Insurance.status + DataSiloName + Neuter_Status + VetCompassBreed, family = "binomial",
                  data = CCL_dat)
pd.cens <- 1 - predict(denom.cens, type = "response")
numer.cens <-
  glm(Lameness_3mo_cens ~ Management.of.CCL.rupture_01, family = binomial(), data = CCL_dat)
pn.cens <- 1 - predict(numer.cens, type = "response")
CCL_dat$sw.c <- pn.cens / pd.cens

#Use combined inverse probability of treatment and censoring weights in binary logistic regression outcome model
library(survey)
#Lameness at approx 3 months
clus <- svydesign(id =~ 1, weights =~ weight$ipw.weights*CCL_dat$sw.c, data = CCL_dat)
res <- svyglm(Lameness_binary_3mo ~ Management.of.CCL.rupture_01 + VetCompassBreed, design = clus)
summary(res)
coef(res)
confint(res)

#Predict Y1 and Y0 (i.e. risk of short-term lameness as if all dogs surgically treated and risk of short-term lameness as if all dogs non-surgically treated)
#Calcuate ATE (Average Treatment Effect)
exp.CCL_dat <- unexp.CCL_dat <- CCL_dat
exp.CCL_dat$Management.of.CCL.rupture_01 <- 1
unexp.CCL_dat$Management.of.CCL.rupture_01 <- 0
Y1 <- predict(res, newdata=exp.CCL_dat, type = "response")
Y0 <- predict(res, newdata=unexp.CCL_dat, type = "response")
ATE.reg <- mean(Y1 -Y0)
ATE.reg;mean(Y1);mean(Y0)

#Bootstrapping the procedure to obtain 95% confidence intervals for the ATE
fATE_marginal <-function(formula, data, indices){
  d <- CCL_dat[indices,]
  exp.CCL_dat <- unexp.CCL_dat <- d
  # set Management.of.CCL.rupture_01=1 in the exposed data and Management.of.CCL.rupture_01=0 in the unexposed data
  exp.CCL_dat$Management.of.CCL.rupture_01 <-1
  unexp.CCL_dat$Management.of.CCL.rupture_01 <- 0
  res <-glm(as.formula(formula), data=d, family="binomial")
  Y1<- predict(res, newdata=exp.CCL_dat, type="response")
  Y0<- predict(res, newdata=unexp.CCL_dat, type="response")
  # our point estimate of risk difference
  ATE_marginal <- mean(Y1 -Y0)
  meanY1=mean(Y1)
  meanY0=mean(Y0)
  return(list=c(ATE_marginal,meanY1,meanY0 ))
  }
#Bootstrapping with 1000 replications
results.marginal <- boot(formula=res, data=CCL_dat, statistic=fATE_marginal,R=1000)
#Get 95% confidence interval
boot.ci.ace.marginal<-boot.ci(results.marginal , conf = 0.95, type = "perc")
boot.ci(results.marginal,type="norm",index=1)

#Calculate risk ratio to use in online E-value calculator (Mathur MB, Ding P, Riddell CA, VanderWeele TJ (2018))
#(using Y1 and Y0 to calculate number of lame dogs exposed etc)
tab <- matrix(c(407,208, 81, 119),byrow=TRUE,nrow=2)
epitab(tab,method="riskratio")

#Long-term lameness outcome model
#Generate IPCW for lameness at 12 months
denom.cens <- glm(Lameness_12mo_cens ~ Bodyweightcat + Bodyweight..kg._missing + Age.at.diagnosis..years. 
                  + Overweight.obese + Overweight.obese_missing + Orthopaedic.comorbidity.at.diagnosis.collapsed + Number.of.non.orthopaedic.comorbidities.at.diagnosis_grouped
                  + Insurance.status + DataSiloName + Neuter_Status + VetCompassBreed, family = "binomial",
                  data = CCL_dat)
pd.cens <- 1 - predict(denom.cens, type = "response")
numer.cens <-
  glm(Lameness_12mo_cens ~ Management.of.CCL.rupture_01, family = binomial(), data = CCL_dat)
pn.cens <- 1 - predict(numer.cens, type = "response")
CCL_dat$sw.c12 <- pn.cens / pd.cens


#Lameness at approx 12 months binary logistic regression outcome model
clus2 <- svydesign(id =~ 1, weights =~ weight$ipw.weights*CCL_dat$sw.c12, data = CCL_dat)
res_2 <- svyglm(Lameness_binary_12mo ~ Management.of.CCL.rupture_01 + VetCompassBreed, design = clus)
summary(res_2)
coef(res_2)
confint(res_2)

#Predict Y1 and Y0
Y1 <- predict(res_2, newdata=exp.CCL_dat, type = "response")
Y0 <- predict(res_2, newdata=unexp.CCL_dat, type = "response")
ATE.reg.2 <- mean(Y1 - Y0)
ATE.reg.2;mean(Y1);mean(Y0)

#Bootstrapping the procedure to obtain 95% confidence intervals for the ATE
fATE_marginal <-function(formula, data, indices){
  d <- CCL_dat[indices,]
  exp.CCL_dat <- unexp.CCL_dat <- d
  # set Management.of.CCL.rupture_01=1 in the exposed data and Management.of.CCL.rupture_01=0 in the unexposed data
  exp.CCL_dat$Management.of.CCL.rupture_01 <-1
  unexp.CCL_dat$Management.of.CCL.rupture_01 <- 0
  res_2 <-glm(as.formula(formula), data=d, family="binomial")
  Y1<- predict(res_2, newdata=exp.CCL_dat, type="response")
  Y0<- predict(res_2, newdata=unexp.CCL_dat, type="response")
  # our point estimate of risk difference
  ATE_marginal <- mean(Y1 -Y0)
  meanY1=mean(Y1)
  meanY0=mean(Y0)
  return(list=c(ATE_marginal,meanY1,meanY0 ))
}

#Bootstrapping with 1000 replications
results.marginal <- boot(formula=res_2, data=CCL_dat, statistic=fATE_marginal,R=1000)
#Get 95% confidence interval
boot.ci.ace.marginal<-boot.ci(results.marginal , conf = 0.95, type = "perc")
boot.ci(results.marginal,type="norm",index=1)

#Calculate risk ratio to use in online E-value calculator
#(using Y1 and Y0 to calculate number of lame dogs exposed etc)
tab <- matrix(c(515,100, 105, 95),byrow=TRUE,nrow=2)
epitab(tab,method="riskratio")

#Analgesia prescription at 3 months outcome model
#Generate IPCW for Analgesia at 3 months
denom.cens <- glm(Analgesia_at_3mo_cens ~ Bodyweightcat + Bodyweight..kg._missing + Age.at.diagnosis..years. 
                  + Overweight.obese + Overweight.obese_missing + Orthopaedic.comorbidity.at.diagnosis.collapsed + Number.of.non.orthopaedic.comorbidities.at.diagnosis_grouped
                  + Insurance.status + DataSiloName + Neuter_Status + VetCompassBreed, family = "binomial",
                  data = CCL_dat)
pd.cens <- 1 - predict(denom.cens, type = "response")
numer.cens <-
  glm(Analgesia_at_3mo_cens ~ Management.of.CCL.rupture_01, family = binomial(), data = CCL_dat)
pn.cens <- 1 - predict(numer.cens, type = "response")
CCL_dat$sw.ca3 <- pn.cens / pd.cens

#Use of analgesia at 3 months binary logistic regression outcome model
clus3 <- svydesign(id =~ 1, weights =~ weight$ipw.weights*CCL_dat$sw.ca3, data = CCL_dat)
res_3 <- svyglm(Analgesia_at_3mo ~ Management.of.CCL.rupture_01 + VetCompassBreed, design = clus)
summary(res_3)
coef(res_3)
confint(res_3)

#Predict Y1 and Y0
Y1 <- predict(res_3, newdata=exp.CCL_dat, type = "response")
Y0 <- predict(res_3, newdata=unexp.CCL_dat, type = "response")
ATE.reg.3 <- mean(Y1 - Y0)
ATE.reg.3;mean(Y1);mean(Y0)

#Bootstrapping the procedure to obtain 95% confidence intervals for the ATE
fATE_marginal <-function(formula, data, indices){
  d <- CCL_dat[indices,]
  exp.CCL_dat <- unexp.CCL_dat <- d
  # set Management.of.CCL.rupture_01=1 in the exposed data and Management.of.CCL.rupture_01=0 in the unexposed data
  exp.CCL_dat$Management.of.CCL.rupture_01 <-1
  unexp.CCL_dat$Management.of.CCL.rupture_01 <- 0
  res_3 <-glm(as.formula(formula), data=d, family="binomial")
  Y1<- predict(res_3, newdata=exp.CCL_dat, type="response")
  Y0<- predict(res_3, newdata=unexp.CCL_dat, type="response")
  # our point estimate of risk difference
  ATE_marginal <- mean(Y1 -Y0)
  meanY1=mean(Y1)
  meanY0=mean(Y0)
  return(list=c(ATE_marginal,meanY1,meanY0 ))
}

#Bootstrapping with 1000 replications
results.marginal <- boot(formula=res_3, data=CCL_dat, statistic=fATE_marginal,R=1000)
#Get 95% confidence interval
boot.ci.ace.marginal<-boot.ci(results.marginal , conf = 0.95, type = "perc")
boot.ci(results.marginal,type="norm",index=1)

#Calculate risk ratio to use in online E-value calculator
#(using Y1 and Y0 to calculate number of dogs prescribed analgesia exposed etc)
tab <- matrix(c(447,168, 68, 132),byrow=TRUE,nrow=2)
epitab(tab,method="riskratio")

#Analgesia prescription at 6 months outcome model
#Generate IPCW for Analgesia at 6 months
denom.cens <- glm(Analgesia_at_6mo_cens ~ Bodyweightcat + Bodyweight..kg._missing + Age.at.diagnosis..years. 
                  + Overweight.obese + Overweight.obese_missing + Orthopaedic.comorbidity.at.diagnosis.collapsed + Number.of.non.orthopaedic.comorbidities.at.diagnosis_grouped
                  + Insurance.status + DataSiloName + Neuter_Status + VetCompassBreed, family = "binomial",
                  data = CCL_dat)
pd.cens <- 1 - predict(denom.cens, type = "response")
numer.cens <-
  glm(Analgesia_at_6mo_cens ~ Management.of.CCL.rupture_01, family = binomial(), data = CCL_dat)
pn.cens <- 1 - predict(numer.cens, type = "response")
CCL_dat$sw.ca6 <- pn.cens / pd.cens

#Use of analgesia at 6 months binary logistic regression outcome model
clus4 <- svydesign(id =~ 1, weights =~ weight$ipw.weights*CCL_dat$sw.ca6, data = CCL_dat)
res_4 <- svyglm(Analgesia_at_6mo ~ Management.of.CCL.rupture_01 + VetCompassBreed, design = clus)
summary(res_4)
coef(res_4)
confint(res_4)

#Predict Y1 and Y0
Y1 <- predict(res_4, newdata=exp.CCL_dat, type = "response")
Y0 <- predict(res_4, newdata=unexp.CCL_dat, type = "response")
ATE.reg.4 <- mean(Y1 - Y0)
ATE.reg.4;mean(Y1);mean(Y0)

#Bootstrapping the procedure to obtain 95% confidence intervals for the ATE
fATE_marginal <-function(formula, data, indices){
  d <- CCL_dat[indices,]
  exp.CCL_dat <- unexp.CCL_dat <- d
  # set Management.of.CCL.rupture_01=1 in the exposed data and Management.of.CCL.rupture_01=0 in the unexposed data
  exp.CCL_dat$Management.of.CCL.rupture_01 <-1
  unexp.CCL_dat$Management.of.CCL.rupture_01 <- 0
  res_4 <-glm(as.formula(formula), data=d, family="binomial")
  Y1<- predict(res_4, newdata=exp.CCL_dat, type="response")
  Y0<- predict(res_4, newdata=unexp.CCL_dat, type="response")
  # our point estimate of risk difference
  ATE_marginal <- mean(Y1 -Y0)
  meanY1=mean(Y1)
  meanY0=mean(Y0)
  return(list=c(ATE_marginal,meanY1,meanY0 ))
}

#Bootstrapping with 1000 replications
results.marginal <- boot(formula=res_4, data=CCL_dat, statistic=fATE_marginal,R=1000)
#Get 95% confidence interval
boot.ci.ace.marginal<-boot.ci(results.marginal , conf = 0.95, type = "perc")
boot.ci(results.marginal,type="norm",index=1)

#Calculate risk ratio to use in online E-value calculator
#(using Y1 and Y0 to calculate number of dogs prescribed analgesia exposed etc)
tab <- matrix(c(486,129, 90, 110),byrow=TRUE,nrow=2)
epitab(tab,method="riskratio")

#Analgesia prescription at 12 months outcome model
#Generate IPCW for Analgesia at 12 months
denom.cens <- glm(Analgesia_at_12mo_cens ~ Bodyweightcat + Bodyweight..kg._missing + Age.at.diagnosis..years. 
                  + Overweight.obese + Overweight.obese_missing + Orthopaedic.comorbidity.at.diagnosis.collapsed + Number.of.non.orthopaedic.comorbidities.at.diagnosis_grouped
                  + Insurance.status + DataSiloName + Neuter_Status + VetCompassBreed, family = "binomial",
                  data = CCL_dat)
pd.cens <- 1 - predict(denom.cens, type = "response")
numer.cens <-
  glm(Analgesia_at_12mo_cens ~ Management.of.CCL.rupture_01, family = binomial(), data = CCL_dat)
pn.cens <- 1 - predict(numer.cens, type = "response")
CCL_dat$sw.ca12 <- pn.cens / pd.cens


#Use of analgesia at 12 months binary logistic regression outcome model
clus5 <- svydesign(id =~ 1, weights =~ weight$ipw.weights*CCL_dat$sw.ca12, data = CCL_dat)
res_5 <- svyglm(Analgesia_at_12mo ~ Management.of.CCL.rupture_01 + VetCompassBreed, design = clus)
summary(res_5)
coef(res_5)
confint(res_5)

#Predict Y1 and Y0
Y1 <- predict(res_5, newdata=exp.CCL_dat, type = "response")
Y0 <- predict(res_5, newdata=unexp.CCL_dat, type = "response")
ATE.reg.5 <- mean(Y1 - Y0)
ATE.reg.5;mean(Y1);mean(Y0)

#Bootstrapping the procedure to obtain 95% confidence intervals for the ATE
fATE_marginal <-function(formula, data, indices){
  d <- CCL_dat[indices,]
  exp.CCL_dat <- unexp.CCL_dat <- d
  # set Management.of.CCL.rupture_01=1 in the exposed data and Management.of.CCL.rupture_01=0 in the unexposed data
  exp.CCL_dat$Management.of.CCL.rupture_01 <-1
  unexp.CCL_dat$Management.of.CCL.rupture_01 <- 0
  res_5 <-glm(as.formula(formula), data=d, family="binomial")
  Y1<- predict(res_5, newdata=exp.CCL_dat, type="response")
  Y0<- predict(res_5, newdata=unexp.CCL_dat, type="response")
  # our point estimate of risk difference
  ATE_marginal <- mean(Y1 -Y0)
  meanY1=mean(Y1)
  meanY0=mean(Y0)
  return(list=c(ATE_marginal,meanY1,meanY0 ))
}

#Bootstrapping with 1000 replications
results.marginal <- boot(formula=res_5, data=CCL_dat, statistic=fATE_marginal,R=1000)
#Get 95% confidence interval
boot.ci.ace.marginal<-boot.ci(results.marginal , conf = 0.95, type = "perc")
boot.ci(results.marginal,type="norm",index=1)

#Calculate risk ratio to use in online E-value calculator
#(using Y1 and Y0 to calculate number of dogs prescribed analgesia exposed etc)
tab <- matrix(c(509,106,100,100),byrow=TRUE,nrow=2)
epitab(tab,method="riskratio")

#Bodyweight as a product term in lameness outcome models to check effect modification
#Lameness at approx 3 months
res_6 <- svyglm(Lameness_binary_3mo ~ Management.of.CCL.rupture*Bodyweightcat + VetCompassBreed, design = clus)
summary(res_6)
coef(res_6)
confint(res_6)

#Lameness at approx 12 months
res_7 <- svyglm(Lameness_binary_12mo ~ Management.of.CCL.rupture*Bodyweightcat + VetCompassBreed, design = clus)
summary(res_7)
coef(res_7)
confint(res_7)

#Calculate Standardised Mean Differences (SMD) to assess balance of covariates
covariateNames <- c("Bodyweightcat", "Age.at.diagnosis..years.", "Overweight.obese",
                    "Orthopaedic.comorbidity.at.diagnosis.collapsed",
                    "Number.of.non.orthopaedic.comorbidities.at.diagnosis_grouped",
                    "Insurance.status", "DataSiloName", "VetCompassBreed", "Neuter_Status")
cclSvy <- svydesign(ids = ~ 1, data = CCL_dat, weights =~ weight$ipw.weights)
tabWeighted <- svyCreateTableOne(vars = covariateNames, strata = "Management.of.CCL.rupture_01", 
                                 data = cclSvy, test = FALSE)
print(tabWeighted, smd = TRUE)

#SMD before weighting
tabUnmatched <- CreateTableOne(vars = covariateNames, strata = "Management.of.CCL.rupture_01", data = CCL_dat, test = FALSE)
## Show table with SMD
print(tabUnmatched, smd = TRUE)

#Summary of IPTW combined with IPCW in different outcome models
summary(weight$ipw.weights*CCL_dat$sw.c)
summary(weight$ipw.weights*CCL_dat$sw.c12)
summary(weight$ipw.weights*CCL_dat$sw.ca3)
summary(weight$ipw.weights*CCL_dat$sw.ca6)
summary(weight$ipw.weights*CCL_dat$sw.ca12)

