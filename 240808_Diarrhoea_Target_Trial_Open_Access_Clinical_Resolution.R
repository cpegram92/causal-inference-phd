# Diarrhoea Target Trial - 21/11/2022
# Research question: "Does antimicrobial prescription compared to no antimicrobial prescription 
# for acute diarrhoea in dogs cause a difference in clinical resolution?"
# Inverse probability of treatment weighting used to adjust for confounding
# Same steps can be applied to gastrointestinal nutraceuticals as the exposure
# Missing indicator approach for missing data

# Packages to use
library(dplyr)
library(forcats)
library(ipw)
library(survey)
library(boot)
library(ggplot2)

# Import raw data (sample data only provided in "230829_Diarrhoea_Target_Trial_Open_Access.csv")
AD_dat <- read.csv("C:/PhD/221130_Diarrhoea_Target_Trial.csv", header = TRUE, stringsAsFactors = TRUE)

# Set factors
factor_vars <- c("Bodyweight_missing", "Insurance_status", "Comorbidity", "Vomiting", "Reduced_appetite", 
                 "Haematochezia", "Pyrexia", "Duration", "Duration_missing", "Veterinary_Group", 
                 "Gastrointestinal_nutraceutical", "Dietary_advice", "Antiparasitic", "Gastrointestinal_agent")

AD_dat[factor_vars] <- lapply(AD_dat[factor_vars], factor)

# Group Breeds
AD_dat$Breed <- fct_other(AD_dat$VetCompass_Breed_Cleaned, 
                          keep = c('Crossbreed', 'Labrador Retriever', 'German Shepherd Dog', 
                                   'Cockapoo', 'Shih-tzu', 'Staffordshire Bull Terrier', 
                                   'French Bulldog', 'Jack Russell Terrier', 'Yorkshire Terrier'),
                          other_level = "Other")

# Categorise continuous variables (-100 = bodyweight not recorded)
AD_dat$Bodyweightcat <- cut(AD_dat$Bodyweight_closest_to_presentation, breaks = c(-101, 0, 10, 20, 30, 90))
AD_dat$Bodyweight_kg <- factor(AD_dat$Bodyweightcat, levels = c("(0,10]", "(10,20]", "(20,30]", "(30,90]", "(-101,0]"))

# Estimate inverse probability of treatment weights
weight <- ipwpoint(
  exposure = Antimicrobial,
  family = "binomial",
  link = "logit",
  numerator = ~ 1,
  denominator = ~ Bodyweight_kg + Bodyweight_missing + Age + Age * Vomiting + I(Age^2) +
    Insurance_status + Comorbidity + Vomiting + Reduced_appetite +
    Haematochezia + Pyrexia + Haematochezia * Pyrexia + Duration + Duration_missing +
    Breed + Veterinary_Group + Gastrointestinal_nutraceutical + Gastrointestinal_nutraceutical * Dietary_advice + 
    Dietary_advice + Antiparasitic + Gastrointestinal_agent,
  data = AD_dat
)

summary(weight$ipw.weights)

# Plot weights
ipwplot(weights = weight$ipw.weights, logscale = FALSE, 
        main = "Stabilized weights", xlim = c(0, 8))

# Using IPW in multivariable logistic regression outcome model
clus <- svydesign(id = ~1, weights = ~weight$ipw.weights, data = AD_dat)
res <- svyglm(Resolution_01 ~ Antimicrobial, design = clus)
summary(res)
coef(res)
confint(res)

# Predict Y1 and Y0 
# (i.e., clinical resolution as if all dogs treated with antimicrobials 
# and clinical resolution as if all dogs not treated with antimicrobials)
# Calcuate ATE (Average Treatment Effect)
exp.AD_dat <- unexp.AD_dat <- AD_dat
exp.AD_dat$Antimicrobial <- 1
unexp.AD_dat$Antimicrobial <- 0
Y1 <- predict(res, newdata = exp.AD_dat, type = "response")
Y0 <- predict(res, newdata = unexp.AD_dat, type = "response")
ATE.reg <- mean(Y1 - Y0)
ATE.reg
mean(Y1)
mean(Y0)

# Bootstrapping the procedure to obtain 95% confidence intervals for the ATE
fATE_marginal <-function(formula, data, indices){
  d <- AD_dat[indices,]
  exp.AD_dat <- unexp.AD_dat <- d
  # set Antimicrobial=1 in the exposed data and Antimicrobial=0 in the unexposed data
  exp.AD_dat$Antimicrobial <-1
  unexp.AD_dat$Antimicrobial <- 0
  res <-glm(as.formula(formula), data=d, family="binomial")
  Y1<- predict(res, newdata=exp.AD_dat, type="response")
  Y0<- predict(res, newdata=unexp.AD_dat, type="response")
  # our point estimate of risk difference
  ATE_marginal <- mean(Y1 -Y0)
  meanY1=mean(Y1)
  meanY0=mean(Y0)
  return(list=c(ATE_marginal,meanY1,meanY0 ))
}

#Bootstrapping with 1000 replications
results.marginal <- boot(formula=res, data=AD_dat, statistic=fATE_marginal,R=1000)

#Get 95% confidence interval
boot.ci.ace.marginal<-boot.ci(results.marginal , conf = 0.95, type = "perc")
boot.ci(results.marginal,type="norm",index=1)

# Calculate Standardised Mean Differences (SMD) to assess balance of covariates
AD_dat$final_weights <- weight$ipw.weights

smds <- tidy_smd(AD_dat,
                 .vars = c("Age", "Antiparasitic", "Bodyweight_kg", "Breed", "Comorbidity", 
                           "Dietary_advice", "Duration", "Gastrointestinal_agent",  
                           "Gastrointestinal_nutraceutical", "Haematochezia", "Insurance_status", 
                           "Pyrexia", "Reduced_appetite", "Veterinary_Group", "Vomiting"),
                 .group = Antimicrobial,
                 .wts = final_weights)

# Plot SMDs
ggplot(data = smds, aes(x = variable, y = abs(smd), group = method, color = method)) +  
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0.1, linetype = "solid", color = "darkgrey", size = 0.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
