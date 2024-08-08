# causal-inference-phd
This repository contains data and R code used for studies as part of my PhD titled "VetCompass eClinical Trials (VETs) – Generating Interventional Evidence from Observational Data"

The file "230829_Diarrhoea_Target_Trial_Open_Access" provides sample data used to answer the research question: "Does antimicrobial prescription compared to no antimicrobial prescription and (separately) gastrointestinal nutraceutical prescription compared to no gastrointestinal nutraceutical prescription for acute diarrhoea in dogs cause a difference in clinical resolution and time to treatment escalation?"
The study used Target Trial Emulation and included dogs diagnosed with acute, uncomplicated diarrhoea between January 1, 2019 and December 31, 2019 within the VetCompass database.
Inverse probability of treatment weighting (IPTW) was used to adjust for confounding.

The files "230913_CCL_Target_Trial_Open_Access" and "230913_CCL_Target_Trial_Open_Access" provides sample data and R code used to answer the research question: 
"Does surgical versus non-surgical management of cranial cruciate ligament rupture in dogs cause different outcomes?"
The study used Target Trial Emulation and included dogs diagnosed with CCL rupture between January 1, 2019 and December 31, 2019 within the VetCompass database
IPTW was used to adjust for confounding, with inverse probability of censoring weighting (IPCW) accounting for censored dogs.

The files "240403_UI_Target_Trial_Open_Access" and "240807_UI_Target_Trial_Open_Access" provides sample data and R code used to answer the research question: 
"What is the causal effect of neutering at 3 to < 7 months versus neutering at ≥ 7 to ≤ 18 months on early-onset urinary incontinence (UI) diagnosis (defined as UI diagnosed at < 8.5 years) in bitches?" 
The study included bitches in the VetCompass database born from January 1, 2010, to December 31, 2012, and neutered between 3 and 18 months old. 
IPTW was used to adjust for confounding, with IPCW accounting for bitches lost to follow-up.

