# TRIAL EMULATION: PIMOBENDAN WITHIN 6 MONTHS OF GRADE IV MURMUR DIAGNOSIS
# CONGESTIVE HEART FAILURE (CHF) OUTCOME WITH COMPETING EVENT DEATH

#Packages
library(survival)
library(boot)
library(ggplot2)

#Read in data
tab<-read.csv("C:/PhD/240423_MVD_Target_Trial_Comp_Ev.csv",sep=",",header=T)
#Include quadratic and interaction term
tab$AgeSquared <- tab$AgeDiagnosis^2
tab$InteractionTerm <- tab$ChronicComorbidityBaseline * tab$AgeDiagnosis

#Clone and censor code dervied from Maringe et al. (2020) tutorial paper

########################################################################
#STEP 1-CLONING: CREATION OF OUTCOME AND FOLLOW-UP TIME IN EACH ARM
########################################################################

#We will create two new variables:
# - fup: the follow-up time in the emulated trial (which can be different from the observed follow-up time)
# - compev: the compev in the emulated trial 

#THESE VARIABLES WILL BE THE OUTCOME AND FOLLOW UP TIME IN THE OUTCOME MODEL

#########################################################################################################################################
#Arm "Control": no treatment within 6 months
tab_control<-tab  # We create a first copy of the dataset: "clones" assigned to the control (no pimobendan) arm
tab_control$arm<-"Control"

#Case 1: Dogs receive Pimobendan within 6 months (of grade IV heart murmur diagnosis): 
#they are still alive and followed-up until pimobendan treatment
tab_control$compev[tab_control$PimoPresc_01==1 & tab_control$TimeDiagnosisPimo <=182.62]<-0

tab_control$fup[tab_control$PimoPresc_01==1
                & tab_control$TimeDiagnosisPimo <=182.62]<-tab_control$TimeDiagnosisPimo[tab_control$PimoPresc_01==1
                                                                                         & tab_control$TimeDiagnosisPimo <=182.62]

#Case 2: Dogs do not receive Pimobendan within 6 months (either no Pimobendan or after 6 months): 
#we keep their observed outcomes and follow-up times 
tab_control$compev[tab_control$PimoPresc_01==0 
                    | (tab_control$PimoPresc_01==1 
                       & tab_control$TimeDiagnosisPimo >182.62)]<-  tab_control$compev[tab_control$PimoPresc_01==0 
                                                                                             | (tab_control$PimoPresc_01==1 
                                                                                                & tab_control$TimeDiagnosisPimo >182.62)]
tab_control$fup[tab_control$PimoPresc_01==0 
                | (tab_control$PimoPresc_01==1 
                   & tab_control$TimeDiagnosisPimo >182.62)]<-  tab_control$FollowUp[tab_control$PimoPresc_01==0 
                                                                                          | (tab_control$PimoPresc_01==1 
                                                                                             & tab_control$TimeDiagnosisPimo >182.62)]

#########################################################################################################################################



#########################################################################################################################################
#Arm "Treated": Pimobendan treatment within 6 months
tab_treated<-tab # We create a second copy of the dataset: "clones" assigned to the treatment arm
tab_treated$arm<-"Treated"

#Case 1: Dogs receive pimobendan within 6 months 
#we keep their observed outcomes and follow-up times
tab_treated$compev[tab_treated$PimoPresc_01==1
                    & tab_treated$TimeDiagnosisPimo <=182.62]<-tab_treated$compev[tab_treated$PimoPresc_01==1
                                                                                        & tab_treated$TimeDiagnosisPimo <=182.62]

tab_treated$fup[tab_treated$PimoPresc_01==1
                & tab_treated$TimeDiagnosisPimo <=182.62]<-tab_treated$FollowUp[tab_treated$PimoPresc_01==1
                                                                                     & tab_treated$TimeDiagnosisPimo <=182.62]

#Case 2: Dogs die or are lost to follow-up before 6 months without having treatment: 
#we keep their observed outcomes and follow-up times
tab_treated$compev[tab_treated$PimoPresc_01==0 
                    & tab_treated$FollowUp <=182.62]<-tab_treated$compev[tab_treated$PimoPresc_01==0 
                                                                                    & tab_treated$FollowUp <=182.62]

tab_treated$fup[tab_treated$PimoPresc_01==0 
                & tab_treated$FollowUp <=182.62]<-tab_treated$FollowUp[tab_treated$PimoPresc_01==0 
                                                                                 & tab_treated$FollowUp <=182.62]

#Case 3: Dogs do not receive treatment within 6 months and are still alive or at risk at 6 months
# they are considered alive and their follow-up time is 6 months
tab_treated$compev[(tab_treated$PimoPresc_01==0 & tab_treated$FollowUp >182.62)
                    | (tab_treated$PimoPresc_01==0  & tab_treated$TimeDiagnosisPimo >182.62)]<-0

tab_treated$fup[(tab_treated$PimoPresc_01==0 & tab_treated$FollowUp >182.62)
                | (tab_treated$PimoPresc_01==0 & tab_treated$TimeDiagnosisPimo >182.62)]<-182.62

#########################################################################################################################################


#################################################################
# STEP1-CLONING: CENSORING STATUS AND FOLLOW-UP TIME UNCENSORED
#################################################################

#We will create two new variables:
# - fup_uncensored: the follow-up time uncensored in the trial arm (can be shorter than the follow-up time in the outcome model)
# - censoring: a binary variable indicating whether the patient was censored in a given arm (either because they receive pimobendan treatment 
#              in the control arm or they didn't receive pimobendan treatment in the treatment arm)

#THESE VARIABLES WILL BE THE OUTCOME AND FOLLOW UP TIME IN THE WEIGHT MODEL

#########################################################################################################################################
#Arm "Control": no treatment within 6 months

#Case 1: Dogs receive pimobendan treatment within 6 months: 
#they are censored in the control group at time of treatment
tab_control$censoring[tab_control$PimoPresc_01==1 & tab_control$TimeDiagnosisPimo <=182.62]<-1

tab_control$fup_uncensored[tab_control$PimoPresc_01==1 
                           & tab_control$TimeDiagnosisPimo <=182.62]<-(tab_control$TimeDiagnosisPimo[tab_control$PimoPresc_01==1 
                                                                                                     & tab_control$TimeDiagnosisPimo <=182.62])

#Case 2: Dogs die or are lost to follow-up before 6 months
#we keep their follow-up time but they are uncensored
tab_control$censoring[tab_control$PimoPresc_01==0 & tab_control$FollowUp <=182.62]<-0

tab_control$fup_uncensored[tab_control$PimoPresc_01==0 & 
                             tab_control$FollowUp <=182.62]<-tab_control$FollowUp[tab_control$PimoPresc_01==0 & 
                                                                                              tab_control$FollowUp <=182.62]

#Case 3: Dogs do not receive pimobendan treatment within 6 months and are still alive or at risk at 6 months
# they are considered uncensored and their follow-up time is 6 months
tab_control$censoring[(tab_control$PimoPresc_01==0 & tab_control$FollowUp >182.62)
                      | (tab_control$PimoPresc_01==1  & tab_control$TimeDiagnosisPimo >182.62)]<-0

tab_control$fup_uncensored[(tab_control$PimoPresc_01==0 & tab_control$FollowUp >182.62)
                           | (tab_control$PimoPresc_01==1  & tab_control$TimeDiagnosisPimo >182.62)]<- 182.62



#########################################################################################################################################


#########################################################################################################################################
#Arm "Treatment": Pimobendan treatment within 6 months

#Case 1: Dogs receive pimobendan treatment within 6 months
# they are uncensored in the treatment arm and remain at risk of censoring until time of treatment
tab_treated$censoring[tab_treated$PimoPresc_01==1 & tab_treated$TimeDiagnosisPimo <=182.62]<-0

tab_treated$fup_uncensored[tab_treated$PimoPresc_01==1 
                           & tab_treated$TimeDiagnosisPimo <=182.62]<-(tab_treated$TimeDiagnosisPimo[tab_treated$PimoPresc_01==1 
                                                                                                     & tab_treated$TimeDiagnosisPimo <=182.62])

#Case 2: Dogs die or are lost to follow-up before 6 months
#we keep their follow-up times but they are uncensored
tab_treated$censoring[tab_treated$PimoPresc_01==0 & tab_treated$FollowUp <=182.62]<-0

tab_treated$fup_uncensored[tab_treated$PimoPresc_01==0 
                           & tab_treated$FollowUp <=182.62]<-tab_treated$FollowUp[tab_treated$PimoPresc_01==0  
                                                                                            & tab_treated$FollowUp <=182.62]

#Case 3: Dogs do not receive treatment within 6 months and are still alive or at risk at 6 months
# they are considered censored and their follow-up time is 6 months
tab_treated$censoring[(tab_treated$PimoPresc_01==0 & tab_treated$FollowUp >182.62)
                      | (tab_treated$PimoPresc_01==1  & tab_treated$TimeDiagnosisPimo >182.62)]<-1

tab_treated$fup_uncensored[(tab_treated$PimoPresc_01==0 & tab_treated$FollowUp >182.62)
                           | (tab_treated$PimoPresc_01==1  & tab_treated$TimeDiagnosisPimo >182.62)]<-182.62


#########################################################################################################################################


###################################################
# CREATION OF THE FINAL DATASET FOR THE ANALYSIS
###################################################

#Each dog appears twice in this dataset (a clone in each treatment arm)
#Combining the two datasets (Control and Treated)
tab<-rbind(tab_control, tab_treated)

tab_old <- tab

# turn outcome to binary for now for the survsplit, change to competing event later on
tab$outcome <- ifelse(tab$compev %in% c(1, 2), 1, tab$compev)

# only select relevant columns
tab <- tab[c("id", "arm","compev","outcome","censoring", "fup","FollowUp", "LTFU","DataSiloName","VetCompassBreed","Insured",
             "ChronicComorbidityBaseline","AgeDiagnosis", "AgeSquared",
             "InteractionTerm","DiagnosticsBaseline")]

# if censoring is 1, make sure LTFU is 0 (artificial censoring is always before LTFU)
tab$LTFU[tab$censoring == 1] <- 0

data_long_Cox <- NULL

##########################
## Bootstrap function
###########################
fboot <- function(tab, indices) {
  total_ids <- unique(tab$id) # get all distinct patients
  total_ids_df <- data.frame(id = total_ids)
  sample_ids <- total_ids_df[indices,] # allows boot to select sample
  
  # select the sample in all arms
  tab <- tab[tab$id %in% sample_ids,] 
  
  ####################################################
  #STEP 2-SPLITTING THE DATASET AT EACH TIME OF EVENT
  ####################################################
  
  
  #Dataframe containing the time of events and an ID for the times of events
  t_events<-sort(unique(tab$fup))
  times<-data.frame("tevent"=t_events,"ID_t"=seq(1:length(t_events)))
  
  
  ####################################
  # Arm "Treated" FIRST
  ####################################
  
  
  tab_t<-tab[tab$arm=="Treated",]
  
  #Creation of the entry variable (Tstart, 0 for everyone)
  tab_t$Tstart<-0
  
  #Splitting the dataset at each time of event until the event happens and sorting it
  data.long<-survSplit(tab_t, cut=t_events, end="fup", 
                       start="Tstart", event="outcome",id="ID") 
  data.long<-data.long[order(data.long$ID,data.long$fup),] 
  
  #Splitting the original dataset at each time of event and sorting it
  #until censoring happens. This is to have the censoring status at each time of event 
  data.long.cens<-survSplit(tab_t, cut=t_events, end="fup", 
                            start="Tstart", event="censoring",id="ID") 
  data.long.cens<-data.long.cens[order(data.long.cens$ID,data.long.cens$fup),] 

  #Splitting the original dataset at each time of event and sorting it
  #until LTFU censoring happens. This is to have the LTFU censoring status at each time of event 
  data.long.cens.LTFU<-survSplit(tab_t, cut=t_events, end="fup", 
                            start="Tstart", event="LTFU",id="ID") 
  data.long.cens.LTFU<-data.long.cens.LTFU[order(data.long.cens.LTFU$ID,data.long.cens.LTFU$fup),] 
  
  #Replacing the censoring variable in data.long by the censoring variable obtained
  # in the second and third split dataset
  data.long$censoring<-data.long.cens$censoring
  data.long$LTFU<-data.long.cens.LTFU$LTFU
  
  #Creating Tstop (end of the interval) 
  data.long$Tstop<-data.long$fup
  
  #Merge and sort
  data.long<-merge(data.long,times,by.x="Tstart",by.y="tevent",all.x=T)
  data.long<-data.long[order(data.long$ID,data.long$fup),] 
  data.long$ID_t[is.na(data.long$ID_t)]<-0
  
  
  ####################################
  # ARM "No Treatment" NOW
  ###################################
  
  
  tab_c<-tab[tab$arm=="Control",]
  
  #Creation of the entry variable (Tstart, 0 for everyone)
  tab_c$Tstart<-0
  
  #Splitting the dataset first at each time of event
  #until the event happens 
  data.long2<-survSplit(tab_c, cut=t_events, end="fup", 
                        start="Tstart", event="outcome",id="ID") 
  data.long2<-data.long2[order(data.long2$ID,data.long2$fup),] 
  
  #Splitting the original dataset at each time of event
  #until censoring happens 
  data.long.cens2<-survSplit(tab_c, cut=t_events, end="fup", 
                             start="Tstart", event="censoring",id="ID") 
  data.long.cens2<-data.long.cens2[order(data.long.cens2$ID,data.long.cens2$fup),]

  #Splitting the original dataset at each time of event and sorting it
  #until LTFU censoring happens. This is to have the LTFU censoring status at each time of event 
  data.long.cens.LTFU2<-survSplit(tab_c, cut=t_events, end="fup", 
                                 start="Tstart", event="LTFU",id="ID") 
  data.long.cens.LTFU2<-data.long.cens.LTFU2[order(data.long.cens.LTFU2$ID,data.long.cens.LTFU2$fup),] 
  
  
  #Replacing the censoring variable in data.long by the censoring variable obtained
  # in the second and third split dataset
  data.long2$censoring<-data.long.cens2$censoring
  data.long2$LTFU<-data.long.cens.LTFU2$LTFU
  
  #Creating Tstop (end of the interval)
  data.long2$Tstop<-data.long2$fup
  
  #Merge and sort
  data.long2<-merge(data.long2,times,by.x="Tstart",by.y="tevent",all.x=T)
  data.long2<-data.long2[order(data.long2$ID,data.long2$fup),] 
  data.long2$ID_t[is.na(data.long2$ID_t)]<-0
  
  #Final dataset
  data<-rbind(data.long,data.long2)
  data_final<-merge(data,times,by="ID_t",all.x=T)
  data_final<-data_final[order(data_final$ID,data_final$fup),]
  
  # if outcome is death (compev=2), change outcome to 2
  data_final$outcome[data_final$outcome == 1 & data_final$compev == 2] <- 2
  
  # only select relevant columns
  data_final <- data_final[c("id", "arm", "ID_t","Tstart","Tstop","outcome", "censoring","LTFU", "fup","DataSiloName","VetCompassBreed","Insured",
                             "ChronicComorbidityBaseline","AgeDiagnosis", "AgeSquared",
                             "InteractionTerm","DiagnosticsBaseline")]
  
  
  ############################################
  #STEP 3- ESTIMATING THE ARTIFICIAL CENSORING WEIGHTS
  ############################################
  
  #######################################################################################################################
  # Arm "Treated" first
  
  data.long<-data_final[data_final$arm=="Treated",]
  
  ###########################
  # STEP 1: censoring model
  ###########################
  
  
  #Cox model
  ms_cens<-coxph(Surv(Tstart, Tstop, censoring)~ DataSiloName + VetCompassBreed +
                   Insured + ChronicComorbidityBaseline + AgeDiagnosis + AgeSquared + 
                   InteractionTerm + DiagnosticsBaseline, ties="efron", data=data.long) 
  
  
  ###########################################################
  # STEP 2: Estimating the probability of remaining uncensored
  ###########################################################
  
  #Design matrix
  design_mat<-as.matrix(data.long[,c("DataSiloName","VetCompassBreed","Insured",
                                     "ChronicComorbidityBaseline","AgeDiagnosis", "AgeSquared",
                                     "InteractionTerm","DiagnosticsBaseline")])
  #Vector of regression coefficients
  beta<-coef(ms_cens)
  
  #Calculation of XB (linear combination of the covariates)
  data.long$lin_pred<-design_mat%*%beta
  
  #Estimating the cumulative hazard (when covariates=0)
  dat.base<-data.frame(basehaz(ms_cens,centered=F))
  names(dat.base)<-c("hazard","t")
  dat.base<-unique(merge(dat.base,times,by.x="t",by.y="tevent",all.x=T))
  
  
  #Merging and reordering the dataset
  data.long<-merge(data.long,dat.base,by="ID_t",all.x=T)
  data.long<-data.long[order(data.long$id,data.long$fup),]
  data.long$hazard<-ifelse(is.na(data.long$hazard),0,data.long$hazard)
  
  
  #Estimating the probability of remaining uncensored at each time of event
  data.long$P_uncens<-exp(-(data.long$hazard)*exp(data.long$lin_pred))  
  
  
  #############################
  # Computing IPC weights
  #############################
  
  #Weights are the inverse of the probability of remaining uncensored
  data.long$weight_Cox<-1/data.long$P_uncens
  ####################################################################################################################
  
  ####################################
  # Arm "No Treatment" now
  
  data.long2<-data_final[data_final$arm=="Control",]
  
  ###########################
  # STEP 1: censoring model
  ###########################
  
  
  #Cox model
  ms_cens2<-coxph(Surv(Tstart, Tstop, censoring) ~ DataSiloName + VetCompassBreed +
                    Insured + ChronicComorbidityBaseline + AgeDiagnosis + AgeSquared + 
                    InteractionTerm + DiagnosticsBaseline, ties="efron", data=data.long2)
  summary(ms_cens2)
  
  
  ###########################################################
  # STEP 2: estimate the probability of remaining uncensored
  ###########################################################
  
  #Design matrix
  design_mat2<-as.matrix(data.long2[,c("DataSiloName","VetCompassBreed","Insured",
                                       "ChronicComorbidityBaseline","AgeDiagnosis","AgeSquared",
                                       "InteractionTerm", "DiagnosticsBaseline")])
  #Vector of regression coefficients
  beta2<-coef(ms_cens2)
  
  #Calculation of XB (linear combineation of the covariates)
  data.long2$lin_pred<-design_mat2%*%beta2
  
  #Estimating the cumulative hazard (when covariates=0)
  dat.base2<-data.frame(basehaz(ms_cens2,centered=F))
  names(dat.base2)<-c("hazard","t")
  
  
  dat.base2<-unique(merge(dat.base2,times,by.x="t",by.y="tevent",all.x=T))
  
  #Merging and reordering the dataset
  data.long2<-merge(data.long2,dat.base2,by="ID_t",all.x=T)
  data.long2<-data.long2[order(data.long2$id,data.long2$fup),]
  data.long2$hazard<-ifelse(is.na(data.long2$hazard),0,data.long2$hazard)
  
  
  #Estimating the probability of remaining uncensored at each time of event
  data.long2$P_uncens<-exp(-(data.long2$hazard)*exp(data.long2$lin_pred))
  
  
  #############################
  # Computing the IPC weights
  #############################
  
  #Weights are the inverse of the probability of remaining uncensored
  data.long2$weight_Cox<-1/data.long2$P_uncens
  data.long2$weight_Cox[data.long2$ID_t==0]<-1
  
  data.long.Cox<-rbind(data.long,data.long2)
  
  # missing outcome means no outcome
  data.long.Cox$outcome[is.na(data.long.Cox$outcome)] <- 0
  

  ############################################
  #STEP 4- ESTIMATING THE LTFU CENSORING WEIGHTS 
  ############################################
  
  #######################################################################################################################
  # Arm "Treated" first
  
  data.long3<-data.long.Cox[data.long.Cox$arm=="Treated",]
  
  ###########################
  # STEP 1: LTFU censoring model
  ###########################
  
  
  #Cox model
  ms_cens3<-coxph(Surv(Tstart, Tstop, LTFU)~ DataSiloName + VetCompassBreed +
                   Insured + ChronicComorbidityBaseline + AgeDiagnosis + AgeSquared + 
                   DiagnosticsBaseline, ties="efron", data=data.long3) 
  
  
  ###########################################################
  # Estimating the probability of remaining uncensored
  ###########################################################
  
  #Design matrix
  design_mat3<-as.matrix(data.long3[,c("DataSiloName","VetCompassBreed","Insured",
                                     "ChronicComorbidityBaseline","AgeDiagnosis", "AgeSquared",
                                     "DiagnosticsBaseline")])
  #Vector of regression coefficients
  beta3<-coef(ms_cens3)
  
  #Calculation of XB (linear combination of the covariates)
  data.long3$lin_pred_LTFU<-design_mat3%*%beta3
  
  #Estimating the cumulative hazard (when covariates=0)
  dat.base3<-data.frame(basehaz(ms_cens3,centered=F))
  names(dat.base3)<-c("hazard_LTFU","t_LTFU")
  dat.base3<-unique(merge(dat.base3,times,by.x="t_LTFU",by.y="tevent",all.x=T))
  
  #Merging and reordering the dataset
  data.long3<-merge(data.long3,dat.base3,by="ID_t",all.x=T)
  data.long3<-data.long3[order(data.long3$id,data.long3$fup),]
  data.long3$hazard_LTFU<-ifelse(is.na(data.long3$hazard_LTFU),0,data.long3$hazard_LTFU)
  
  #Estimating the probability of remaining uncensored at each time of event
  data.long3$P_uncens_LTFU<-exp(-(data.long3$hazard_LTFU)*exp(data.long3$lin_pred_LTFU))  
  
  
  #############################
  # Computing IPC weights
  #############################
  
  #Weights are the inverse of the probability of remaining uncensored
  data.long3$weight_Cox_LTFU<-1/data.long3$P_uncens_LTFU
  ####################################################################################################################
  
  ####################################
  # Arm "No Treatment" now
  
  data.long4<-data.long.Cox[data.long.Cox$arm=="Control",]
  
  ###########################
  # STEP 1: censoring model
  ###########################
  
  
  #Cox model
  ms_cens4<-coxph(Surv(Tstart, Tstop, LTFU) ~ DataSiloName + VetCompassBreed +
                    Insured + ChronicComorbidityBaseline + AgeDiagnosis + AgeSquared + 
                    DiagnosticsBaseline, ties="efron", data=data.long4)
  summary(ms_cens4)
  
  
  ###########################################################
  # STEP 2: estimate the probability of remaining uncensored
  ###########################################################
  
  #Design matrix
  design_mat4<-as.matrix(data.long4[,c("DataSiloName","VetCompassBreed","Insured",
                                       "ChronicComorbidityBaseline","AgeDiagnosis","AgeSquared",
                                       "DiagnosticsBaseline")])
  #Vector of regression coefficients
  beta4<-coef(ms_cens4)
  
  #Calculation of XB (linear combineation of the covariates)
  data.long4$lin_pred_LTFU<-design_mat4%*%beta4
  
  #Estimating the cumulative hazard (when covariates=0)
  dat.base4<-data.frame(basehaz(ms_cens4,centered=F))
  names(dat.base4)<-c("hazard_LTFU","t_LTFU")
  
  
  dat.base4<-unique(merge(dat.base4,times,by.x="t_LTFU",by.y="tevent",all.x=T))
  
  #Merging and reordering the dataset
  data.long4<-merge(data.long4,dat.base4,by="ID_t",all.x=T)
  data.long4<-data.long4[order(data.long4$id,data.long4$fup),]
  data.long4$hazard_LTFU<-ifelse(is.na(data.long4$hazard_LTFU),0,data.long4$hazard_LTFU)
  
  
  #Estimating the probability of remaining uncensored at each time of event
  data.long4$P_uncens_LTFU<-exp(-(data.long4$hazard_LTFU)*exp(data.long4$lin_pred_LTFU))
  
  
  #############################
  # Computing the IPC weights
  #############################
  
  
  #Weights are the inverse of the probability of remaining uncensored
  data.long4$weight_Cox_LTFU<-1/data.long4$P_uncens_LTFU
  data.long4$weight_Cox_LTFU[data.long4$ID_t==0]<-1
  
  data.long.Cox2<-rbind(data.long3,data.long4)
  
  # missing outcome means no outcome
  data.long.Cox2$outcome[is.na(data.long.Cox2$outcome)] <- 0
  
  
  ############################################
  # Obtain one final weight: artificial censoring weight * LTFU weight
  ############################################
  # Calculate the final weight
  data.long.Cox2$weight_final <- data.long.Cox2$weight_Cox * data.long.Cox2$weight_Cox_LTFU
  # Calculate the 1st and 99th percentiles of weight_final
  percentile_1 <- quantile(data.long.Cox2$weight_final, 0.01)
  percentile_99 <- quantile(data.long.Cox2$weight_final, 0.99)
  
  # Apply truncation at percentiles
  data.long.Cox2$weight_final <- pmin(pmax(data.long.Cox2$weight_final, percentile_1), percentile_99)
  
  ############################################
  # STEP 5- ESTIMATING THE SURVIVOR FUNCTION
  ############################################
  
  
  ###################################################
  # Aalen Johansen
  ###################################################
  
  data.long.Cox2$outcome <- factor(data.long.Cox2$outcome, 0:2, c("atrisk", "CHF", "death"))
  
  ci_Treated <- survfit(Surv(Tstart, Tstop, outcome) ~ 1, data = data.long.Cox2[data.long.Cox2$arm=='Treated',], weights = weight_final[,1], id = id)
  ci_Control <- survfit(Surv(Tstart, Tstop, outcome) ~ 1, data = data.long.Cox2[data.long.Cox2$arm=='Control',], weights = weight_final[,1], id = id)
  
  time_points_plot <- sort(unique(data.long.Cox2$Tstart))
  
  # outcome of CHF, accounting for competing event death
  ci_Treated_CHF <- summary(ci_Treated, times = time_points_plot, extend=TRUE)$pstate[,2]
  ci_Control_CHF <-summary(ci_Control, times = time_points_plot, extend=TRUE)$pstate[,2]
  
  df_CHF <- data.frame(time_points = time_points_plot, treated = ci_Treated_CHF, control = ci_Control_CHF)
  
  # competing event death 
  ci_Treated_death <- summary(ci_Treated, times = time_points_plot, extend=TRUE)$pstate[,3]
  ci_Control_death <-summary(ci_Control, times = time_points_plot, extend=TRUE)$pstate[,3]
  
  df_death <- data.frame(time_points = time_points_plot, treated = ci_Treated_death, control = ci_Control_death)
  
  # Assign the plot dataframes to a variable in the global environment
  assign("df_CHF", df_CHF, envir = .GlobalEnv)
  assign("df_death", df_death, envir = .GlobalEnv)
  
  # Calculate 5-year RMTL and difference in RMTL
  rmtl_Treated <- summary(ci_Treated, rmean=1826)$table[2, 3]
  rmtl_Control <- summary(ci_Control, rmean=1826)$table[2, 3]
  RMTL_diff <- rmtl_Treated - rmtl_Control
  
  # find 5-year risk of CHF
  ci_Treated_risk_CHF <- summary(ci_Treated, times=1826, extend=TRUE)$pstate[,2]
  ci_Control_risk_CHF <- summary(ci_Control, times=1826, extend=TRUE)$pstate[,2]
  
  # find 5-year risk of death
  ci_Treated_risk_death <- summary(ci_Treated, times=1826, extend=TRUE)$pstate[,3]
  ci_Control_risk_death <- summary(ci_Control, times=1826, extend=TRUE)$pstate[,3]
  
  # risk difference for CHF and death
  risk_diff_CHF <- ci_Treated_risk_CHF - ci_Control_risk_CHF
  risk_diff_death <- ci_Treated_risk_death - ci_Control_risk_death
  
  res <- c(ci_Treated_risk_CHF, ci_Control_risk_CHF, risk_diff_CHF, ci_Treated_risk_death, ci_Control_risk_death, risk_diff_death, rmtl_Treated, rmtl_Control, RMTL_diff)
  names(res) <- c("5_year_risk_CHF_treated", "5_year_risk_CHF_control", "5_year_risk_difference_CHF", "5_year_risk_death_treated", "5_year_risk_death_control", "5_year_risk_difference_death", "RMTL_Treated", "RMTL_Control", "RMTL_difference")
  
  return(res)
}

##############################################
# CALLING THE BOOTSTRAP FUNCTION
##############################################

# Bootstrapping with 100 replications
results <- boot(data=tab, statistic=fboot, R=1000)
results

# 95% confidence intervals for each measure
boot.ci(results, type="norm", index=1) # 5_year_risk_CHF_Treated
boot.ci(results, type="norm", index=2) # 5_year_risk_CHF_control
boot.ci(results, type="norm", index=3) # 5_year_risk_difference_CHF
boot.ci(results, type="norm", index=4) # 5_year_risk_death_treated
boot.ci(results, type="norm", index=5) # 5_year_risk_death_control
boot.ci(results, type="norm", index=6) # 5_year_risk_difference_death
boot.ci(results, type="norm", index=7) # RMTL_Treated
boot.ci(results, type="norm", index=8) # RMTL_Control
boot.ci(results, type="norm", index=9) # RMTL_difference

##############################################
# PLOTTING CUMULATIVE INCIDENCE CURVES
##############################################

# Cumulative incidence curves based on the Aalen-Johanson estimator
# CHF
ggplot(df_CHF, aes(x = time_points)) +
  geom_line(aes(y = treated, color = "Pimobendan Prescribed"), size = 1.3) +#, color = "darkgreen"
  geom_line(aes(y = control, color = "Pimobendan Not Prescribed"), size = 1.3) +#, color = "deeppink2"
  labs(x = "Days since grade IV murmur diagnosis", y = "Cumulative incidence of CHF", color = "") +
  theme_minimal() +
  scale_color_manual(values = c("Pimobendan Prescribed"="darkgreen", "Pimobendan Not Prescribed"="deeppink2")) +
  theme(legend.position = "top")

# competing event death
ggplot(df_death, aes(x = time_points)) +
  geom_line(aes(y = treated, color = "Pimobendan Prescribed"), size = 1.3) +#, color = "darkgreen"
  geom_line(aes(y = control, color = "Pimobendan Not Prescribed"), size = 1.3) +#, color = "deeppink2"
  labs(x = "Days since grade IV murmur diagnosis", y = "Cumulative incidence of competing event death", color = "") +
  theme_minimal() +
  scale_color_manual(values = c("Pimobendan Prescribed"="darkgreen", "Pimobendan Not Prescribed"="deeppink2")) +
  theme(legend.position = "top")

# Combine the two data frames
df_combined <- merge(df_CHF, df_death, by = "time_points", suffixes = c("_CHF", "_death"))

# Combine the two data frames
df_combined <- merge(df_CHF, df_death, by = "time_points", suffixes = c("_CHF", "_death"))

# Plot combined graph
ggplot(df_combined, aes(x = time_points)) +
  geom_line(aes(y = treated_CHF, color = "CHF - Pimobendan Prescribed"), size = 1.3, linetype = "solid") +
  geom_line(aes(y = control_CHF, color = "CHF - Pimobendan Not Prescribed"), size = 1.3, linetype = "solid") +
  geom_line(aes(y = treated_death, color = "Death - Pimobendan Prescribed"), size = 0.8, linetype = "dashed") +
  geom_line(aes(y = control_death, color = "Death - Pimobendan Not Prescribed"), size = 0.8, linetype = "dashed") +
  labs(x = "Days since grade IV murmur diagnosis", y = "Cumulative incidence", color = "") +
  scale_color_manual(values = c("CHF - Pimobendan Prescribed" = "darkgreen",
                                "CHF - Pimobendan Not Prescribed" = "deeppink2",
                                "Death - Pimobendan Prescribed" = "darkgreen",
                                "Death - Pimobendan Not Prescribed" = "deeppink2")) +
  theme_minimal() +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1))


