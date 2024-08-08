# UI Target Trial Sample R Code - 16/10/2023
# Causal effect of neutering at 3 to < 7 months versus neutering at ≥ 7 to ≤ 18 months 
# on early-onset urinary incontinence (UI) diagnosis in bitches.
# Inverse probability of treatment weighting (IPTW) used to adjust for confounding 
# and inverse probability of censoring weighting (IPCW) used to address loss to follow-up.

# Load required packages
library(dplyr)
library(ggplot2)
library(forcats)
library(ipw)
library(boot)
library(clubSandwich)
library(sandwich)
library(finalfit)

# Import raw data (sample data only provided in "240610_UI_Target_Trial_Open_Access.csv")
UI_dat <- read.csv("C:/PhD/240219_UI_Cases_And_Non_Cases.csv", header = TRUE, stringsAsFactors = TRUE)

# Set factors
factor_vars <- c("Veterinary_Group", "Insurance_status", "Comorbidity", "Umbilical_hernia", 
                 "Skin_disorder", "Ear_disorder", "Urogenital_disorder", "UI", 
                 "Neuter_Later_01", "Time_interval", "AtRiskBreed")

UI_dat[factor_vars] <- lapply(UI_dat[factor_vars], factor)

# Group breeds
UI_dat$Breed <- fct_other(UI_dat$Breed, 
                          keep = c('Crossbreed', 'Labrador Retriever', 'German Shepherd Dog', 
                                   'Shih-tzu', 'Staffordshire Bull Terrier', 'English Springer Spaniel',
                                   'Boxer', 'Border Collie', 'Cocker Spaniel', 'West Highland White Terrier',
                                   'Jack Russell Terrier', 'Yorkshire Terrier'), 
                          other_level = "Other")

# Calculate inverse probability of treatment weights (IPTW)
weight <- ipwpoint(
  exposure = Neuter_Later_01,
  family = "binomial",
  link = "logit",
  numerator = ~ 1,
  denominator = ~ Breed + Veterinary_Group + Insurance_status + Comorbidity + 
    Umbilical_hernia + Skin_disorder + Ear_disorder + Urogenital_disorder + 
    Breed * Comorbidity,
  data = UI_dat
)

summary(weight$ipw.weights)

# Plot weights
ipwplot(weights = weight$ipw.weights, logscale = FALSE, 
        main = "Stabilized weights", xlim = c(0, 8))

# Cumulative inverse probability of censoring weights (Cens_wt) calculated in separate R code and imported via CSV
# Truncate weights by capping within the 1st and 99th percentiles
percentiles <- quantile(UI_dat$Cens_wt, probs = c(0.01, 0.99))
UI_dat$cens_weights_truncated <- pmin(pmax(UI_dat$Cens_wt, percentiles[1]), percentiles[2])

# Summary of overall truncated IPCWs
summary(UI_dat$cens_weights_truncated)

# Multiply IPTW by truncated IPCW to obtain final weights
UI_dat$final_weights <- weight$ipw.weights * UI_dat$cens_weights_truncated

# Summary of overall weights
summary(UI_dat$final_weights)

# Calculate Standardised Mean Differences (SMD) to assess balance of covariates
smds <- tidy_smd(UI_dat,
                 .vars = c("Breed", "Veterinary_Group", "Insurance_status", "Comorbidity", 
                           "Umbilical_hernia", "Skin_disorder", "Ear_disorder", "Urogenital_disorder"),
                 .group = Neuter_Later_01,
                 .wts = final_weights
)

# Plot SMDs
Fig_2 <- ggplot(data = smds, aes(x = variable, y = abs(smd), group = method, color = method)) +  
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0.1, linetype = "solid", color = "darkgrey", size = 0.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate x-axis text

# Fit pooled logistic regression model including IPTW and IPCW
model_ipw <- glm(UI ~ Neuter_Later_01 + Age_at_neuter_months, 
                 data = UI_dat, family = binomial(link = "logit"), 
                 weights = final_weights, na.action = "na.exclude")

summary(model_ipw)

# Calculate odds ratios
odds_ratio <- exp(coef(model_ipw))
print(odds_ratio)

# Remove rows with missing values for the "UI" variable and re-run analysis with cluster-robust variance estimation
UI_dat <- UI_dat[complete.cases(UI_dat$UI), ]

model_ipw <- glm(UI ~ Neuter_Later_01 + Age_at_neuter_months, 
                 data = UI_dat, family = binomial(link = "logit"), 
                 weights = final_weights, na.action = "na.exclude")

vcovCR(model_ipw, cluster = UI_dat$PatientID, type = "CR2")

# Function for bootstrapped 95% confidence intervals with cluster-robust variance
Causal_odds <- function(data, indices) {
  d <- data[indices, ]
  model_ipw <- glm(UI ~ Neuter_Later_01 + Age_at_neuter_months, 
                   data = d, family = binomial(link = "logit"), 
                   weights = d$final_weights, na.action = "na.exclude")
  
  coef_neuter <- coef(model_ipw)[2]
  vcov_cluster <- vcovCR(model_ipw, cluster = d$PatientID, type = "CR2")
  se <- sqrt(vcov_cluster[2, 2])
  odds_ratio <- exp(coef_neuter)
  
  return(odds_ratio)
}

# Bootstrapping with 1000 replications
results.odds <- boot(data = UI_dat, statistic = Causal_odds, R = 1000)

# Get 95% confidence intervals
boot.ci.odds <- boot.ci(results.odds, conf = 0.95, type = "basic")

# Print results
print(results.odds)
print(boot.ci.odds)
