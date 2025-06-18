# causal-inference-phd
This repository contains sample data and R code used for studies as part of my PhD titled "VetCompass eClinical Trials (VETs) – Generating Interventional Evidence from Observational Data"

The files "230829_Diarrhoea_Target_Trial_Open_Access", "240808_Diarrhoea_Target_Trial_Open_Access_Clinical_Resolution" and "240807_Diarrhoea_Target_Trial_Open_Access_Time_To_Event" 
provide sample data and R code used to answer the research question: "Does antimicrobial prescription compared to no antimicrobial prescription and (separately) gastrointestinal nutraceutical prescription 
compared to no gastrointestinal nutraceutical prescription for acute diarrhoea in dogs cause a difference in clinical resolution and time to treatment escalation?"
The study used Target Trial Emulation and included dogs diagnosed with acute, uncomplicated diarrhoea between January 1, 2019 and December 31, 2019 within the VetCompass database.
Inverse probability of treatment weighting (IPTW) was used to adjust for confounding.
The R files show example code for antimicrobial prescription as the exposure.

Published paper reference: Target Trial Emulation: Do antimicrobials or gastrointestinal nutraceuticals prescribed at first presentation for acute diarrhoea cause a better clinical outcome in dogs under primary veterinary care in the UK? Pegram C., Diaz-Ordaz K., Brodbelt D. C., Chang Y., Tayler S., Allerton F., Prisk L., Church D. B. & O'Neill D. G. (2023) PLOS ONE

The files "230913_CCL_Target_Trial_Open_Access" and "230913_CCL_Target_Trial_Open_Access" provides sample data and R code used to answer the research question: 
"Does surgical versus non-surgical management of cranial cruciate ligament rupture in dogs cause different outcomes?"
The study used Target Trial Emulation and included dogs diagnosed with CCL rupture between January 1, 2019 and December 31, 2019 within the VetCompass database
IPTW was used to adjust for confounding, with inverse probability of censoring weighting (IPCW) accounting for censored dogs.

Published paper reference: Target Trial Emulation: Does surgical versus non-surgical management of cranial cruciate ligament rupture in dogs cause different outcomes? Pegram C., Diaz-Ordaz K., Brodbelt, D. C., Chang Y.,  Frykfors von Hekkel A., Wu C., Church, D. B., O’Neill, D. G. (2024) Preventative Veterinary Medicine

The files "240403_UI_Target_Trial_Open_Access" and "240807_UI_Target_Trial_Open_Access" provides sample data and R code used to answer the research question: 
"What is the causal effect of neutering at 3 to < 7 months versus neutering at ≥ 7 to ≤ 18 months on early-onset urinary incontinence (UI) diagnosis (defined as UI diagnosed at < 8.5 years) in bitches?" 
The study included bitches in the VetCompass database born from January 1, 2010, to December 31, 2012, and neutered between 3 and 18 months old. 
IPTW was used to adjust for confounding, with IPCW accounting for bitches lost to follow-up.

Published paper reference: Later-Age Neutering Causes Lower Risk Of Early‐Onset Urinary Incontinence Than Early Neutering – A Vetcompass Target Trial Emulation Study. Pegram, C., Diaz-Ordaz, K., Brodbelt, D. C., Chang, Y., Hall, J. L., Church, D. B. & O’neill D.G. 2024. PLOS ONE.

The files "240807_MVD_Competing_Event_Open_Access" and "250618_MVD_Target_Trial_Comp_Ev_Open_Access" provides sample data and R code used to answer the research question: 
"Does the prescription of pimobendan to dogs within 6 months of first grade IV/VI heart murmur cause a difference in (i) time to congestive heart failure (CHF), considering death as a competing event (primary outcome)?"
The study used Target Trial Emulation and included dogs aged 6 years or older with a bodyweight ≤15 kg, first grade IV/VI heart murmur recorded from January 1st 2016 to December 31st 2018. 
Intention to treat analysis with ‘clone-and-censor’ weighting is used to account for baseline confounding and immortal time bias. The artificial censoring is then addressed using inverse probability weighting. Covariate balance at the end of the grace period is assessed using standardised mean differences. Loss to follow-up beyond the grace period addressed using inverse probability of censoring weighting. A weighted non-parametric estimation of the cumulative probabilities and RMTL using the Aalen-Johansen approach.  Dogs are followed up until artificially censored, CHF, death, loss to follow-up, or end of follow-up period, whichever came first.

Published paper reference: Pegram, C., Diaz-Ordaz, K., Brodbelt, D.C., Chang, Y., Boswood, A., Wilshaw, J., Reep, C.A.T., Balling, S., Sahota, J., Church, D.B., O’Neill D.G. (2025). Emulating the EPIC trial using VetCompass primary-care data: causal effects of pimobendan in UK dogs with grade IV/VI heart murmurs. PLOS One.
