# Diarrhoea Target Trial - 06/12/2022
# Research question: "Does antimicrobial prescription compared to no antimicrobial prescription 
# for acute diarrhoea in dogs cause a difference in time to treatment escalation?"
# Inverse probability of treatment weighting used to adjust for confounding
# Same steps can be applied to gastrointestinal nutraceuticals as the exposure
# Missing indicator approach for missing data

# Load required packages
library(tidyr)
library(dplyr)
library(ggplot2)
library(forcats)
library(ipw)
library(survival)
library(splitstackshape)
library(ggfortify)
library(survminer)

# Import raw data (sample data only provided in "230829_Diarrhoea_Target_Trial_Open_Access.csv")
AD_dat <- read.csv("C:/PhD/221130_Diarrhoea_Target_Trial_Survival.csv", header = TRUE, stringsAsFactors = TRUE)

# Set categorical variables as factors
factors <- c("Bodyweight_missing", "Insurance_Status_01", "Comorbidity.at.diagnosis_01",
             "Vomiting_01", "Reduced.Appetite_01", "Hematochezia_01", "Pyrexia_01",
             "Duration_01", "Duration_missing", "DataSiloName_01", "GN", "Dietary.Management",
             "Antiparasitic", "Gastrointestinal_agent")
AD_dat[factors] <- lapply(AD_dat[factors], factor)

# Group less common breeds into 'Other' category
AD_dat$VetCompassBreed <- fct_other(
  AD_dat$VetCompass_Breed_Cleaned,
  keep = c('Crossbreed', 'Labrador Retriever', 'German Shepherd Dog', 'Cockapoo',
           'Shih-tzu', 'Staffordshire Bull Terrier', 'French Bulldog',
           'Jack Russell Terrier', 'Yorkshire Terrier'),
  other_level = "Other"
)

# Categorise continuous variables
# Bodyweight: -100 indicates bodyweight not recorded
AD_dat$Bodyweightcat <- cut(AD_dat$Bodyweight_closest_to_presentation, breaks = c(-101, 0, 10, 20, 30, 90))
AD_dat$Bodyweight_cat <- factor(AD_dat$Bodyweightcat, levels = c("(0,10]", "(10,20]", "(20,30]", "(30,90]", "(-101,0]"))

# Estimate Inverse Probability of Treatment Weights (IPTW)
weight <- ipwpoint(
  exposure = A,
  family = "binomial",
  link = "logit",
  numerator = ~ 1,
  denominator = ~ Bodyweight_cat + Bodyweight_missing + Age_at_diagnosis + Age_at_diagnosis * Vomiting_01 +
    I(Age_at_diagnosis^2) + Insurance_Status_01 + Comorbidity.at.diagnosis_01 + Vomiting_01 +
    Reduced.Appetite_01 + Hematochezia_01 + Pyrexia_01 + Hematochezia_01 * Pyrexia_01 +
    Duration_01 + Duration_missing + VetCompassBreed + DataSiloName_01 + GN +
    GN * Dietary.Management + Dietary.Management + Antiparasitic + Gastrointestinal_agent,
  data = AD_dat
)

# Summarize weights
summary(weight$ipw.weights)

# Plot the distribution of weights
ipwplot(weights = weight$ipw.weights, logscale = FALSE, main = "Stabilized weights", xlim = c(0, 8))

# Following sections derived from R code by by Joy Shi and Sean McGrath: 
# (https://remlapmot.github.io/cibookex-r/causal-survival-analysis.html)

# Nonparametric estimation of survival curves
fit_2 <- survfit(Surv(AD_dat$Time_to_event, AD_dat$Treatment_escalation) ~ AD_dat$A, data = AD_dat)
ggsurvplot(fit_2, data = AD_dat, xlab = "Days of follow-up", ylab = "Time-to-event probability",
           main = "Product-Limit Survival Estimates", risk.table = TRUE)

# Estimation of survival curves via IP weighted hazards model
# Code excluding 95% confidence intervals for simplicity 

# Estimation of denominator of IP weights using logistic regression
p.denom <- glm(A ~ Bodyweight_cat + Bodyweight_missing + Age_at_diagnosis + Age_at_diagnosis * Vomiting_01 +
                 I(Age_at_diagnosis^2) + Insurance_Status_01 + Comorbidity.at.diagnosis_01 + Vomiting_01 +
                 Reduced.Appetite_01 + Hematochezia_01 + Pyrexia_01 + Hematochezia_01 * Pyrexia_01 +
                 Duration_01 + Duration_missing + VetCompassBreed + DataSiloName_01 + GN +
                 GN * Dietary.Management + Dietary.Management + Antiparasitic + Gastrointestinal_agent, 
               data = AD_dat, family = binomial())
AD_dat$pd.A <- predict(p.denom, AD_dat, type = "response")

# Estimation of numerator of IP weights using logistic regression with intercept only
p.num <- glm(A ~ 1, data = AD_dat, family = binomial())
AD_dat$pn.A <- predict(p.num, AD_dat, type = "response")

# Computation of estimated weights
AD_dat$sw.a <- ifelse(AD_dat$A == 1, AD_dat$pn.A / AD_dat$pd.A, (1 - AD_dat$pn.A) / (1 - AD_dat$pd.A))

# Creation of dog-day data (expanding each row based on Time_to_event)
AD_dat.ipw <- expandRows(AD_dat, "Time_to_event", drop = FALSE) 
AD_dat.ipw$time <- sequence(rle(AD_dat.ipw$PatientID)$lengths) - 1
AD_dat.ipw$event <- ifelse(AD_dat.ipw$time == AD_dat.ipw$Time_to_event - 1 & AD_dat.ipw$Treatment_escalation == 1, 1, 0)
AD_dat.ipw$timesq <- AD_dat.ipw$time^2

# Fit the weighted hazards model
ipw.model <- glm(event == 0 ~ A + I(A * time) + I(A * timesq) + time + timesq, family = binomial(), weight = sw.a, data = AD_dat.ipw)
summary(ipw.model)

# Create data for survival curves
ipw.A_0 <- data.frame(cbind(seq(0, 30), 0, (seq(0, 30))^2))
ipw.A_1 <- data.frame(cbind(seq(0, 30), 1, (seq(0, 30))^2))

colnames(ipw.A_0) <- c("time", "A", "timesq")
colnames(ipw.A_1) <- c("time", "A", "timesq")

# Estimate (1-hazard) for each dog-day
ipw.A_0$p.noevent0 <- predict(ipw.model, ipw.A_0, type = "response")
ipw.A_1$p.noevent1 <- predict(ipw.model, ipw.A_1, type = "response")
ipw.A_0$p.noevent0[1] <- 1  # Set initial survival probability to 1
ipw.A_1$p.noevent1[1] <- 1  # Set initial survival probability to 1

# Compute survival probabilities for each dog-day
ipw.A_0$surv0 <- cumprod(ipw.A_0$p.noevent0)
ipw.A_1$surv1 <- cumprod(ipw.A_1$p.noevent1)

# Prepare data for plotting survival curves
ipw.graph <- merge(ipw.A_0, ipw.A_1, by = c("time", "timesq"))
ipw.graph$survdiff <- ipw.graph$surv1 - ipw.graph$surv0

# Plot the estimated survival curves
ggplot(ipw.graph, aes(x = time, y = surv)) + 
  geom_line(aes(y = surv0, colour = "0")) + 
  geom_line(aes(y = surv1, colour = "1")) + 
  xlab("Days") + 
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 2)) +
  scale_y_continuous(limits = c(0.8, 1), breaks = seq(0.8, 1, 0.2)) +
  ylab("Survival") + 
  ggtitle("Survival from IP weighted hazards model") + 
  labs(colour = "A:") +
  theme_bw() + 
  theme(legend.position = "bottom")
