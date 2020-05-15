###############################################################################################
##                              Analyses relative to the paper:                              ##
##         Interictal spike during AED's withdrawal: indicator of the underlying pathology?  ##
###############################################################################################

#### Authors: 
# Pia De Stefano* (corresponding author)
# Doctor resident in Neurology
# EEG and Epilepsy Unit, University  Hospitals of Geneva.
# Mail : pia.destefano@hcuge.ch
# Tel : 0041(0)795533812
# ORCID: https://orcid.org/0000-0002-7979-0994
# 
# Eric MENETRE
# Research Assistant
# EEG and Epilepsy Unit, University Hospitals of Geneva.
# Mail:Eric.M?n?tr?@hcuge.ch
# ORCID: https://orcid.org/0000-0002-1101-1288
# 
# Pieter van Mierlo
# Medical Image and Signal Processing Group, Department of Electronics and Information Systems, Ghent University, Ghent, Belgium  
# Mail: Pieter.vanMierlo@ugent.be
# 
# Margitta Seeck
# Doctor and Professor, head of the EEG and Epilepsy Unit
# EEG and Epilepsy Unit, University Hospitals of Geneva.
# Mail : margitta.seeck@hcuge.ch
# Tel : 0041(0)795533820
# ORCID: https://orcid.org/0000-0002-6702-0167

##----------------------------------------------------------------
##                            Packages                           -
##----------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(bannerCommenter)
library(MASS)
library(lme4)
library(lmerTest)
library(NPL) # from https://github.com/EricMenetre/NPL -> Follow instruction on the README regarding installation
library(readxl)
library(optimx)
library(ggpubr)
library(MuMIn)
library(car)
library(emmeans)

##---------------------------------------------------------------
##    Changes in spiking activity during the pre-ictal period   -
##---------------------------------------------------------------

# Data import
data_N_sz <- read_excel("data_N_sz.xlsx")
View(data_N_sz)

# Transformation of the data to satisfy the expectations of the tidy format
data_tidy <- data_N_sz%>%
  filter(N_crise_P <= 4)%>%
  dplyr::select(Spikes_Baseline, Spikes_Seizure, Patient, N_crise_P, Type, Localisation, Les)%>%
  gather(key = "time", value = "spikes", -Patient, -N_crise_P, -Type, -Localisation, -Les)%>%
  mutate(spikes_w0 = spikes + 0.0001, # Laplacian transformation of the data
         log.spikes = log(spikes_w0),
         localisation_temp_ext = case_when(Localisation == "R temporal" ~ "temporal",
                                           Localisation == "L temporal" ~ "temporal",
                                           TRUE ~ "extratemporal"))

View(data_tidy)

# Transformation of all the categorical variables to factor
data_tidy$Patient <- factor(data_tidy$Patient)
data_tidy$Type <- factor(data_tidy$Type)
data_tidy$Localisation <- factor(data_tidy$Localisation)
data_tidy$Les <- factor(data_tidy$Les)
data_tidy$time <- factor(data_tidy$time)

str(data_tidy)

# Exploration of the distribution of the DV (number of spikes) and exploration of the variability of the random factors of the model
exp_plots_LMM(data_tidy, data_tidy$spikes, data_tidy$Patient, data_tidy$N_crise_P) # The DV is not normally distributed (more of a Poisson distribution)

# STATISTICAL MODELING
# From the simpler (empty) model to the most complete one, comparaison of which one reduces the most the AIC, BIC and deviance
m0 <- glmer(spikes ~ 1 + (1|Patient),
            family = poisson(link = "log"),
            data = data_tidy) 

summary(m0)
LMM_check(m0)

m1 <- glmer(spikes ~ time + (1|Patient),
            family = poisson(link = "log"),
            data = data_tidy) 

summary(m1)
LMM_check(m1)
anova(m0, m1) # Significant

m2 <- glmer(spikes ~ time + N_crise_P + (1|Patient),
            family = poisson(link = "log"),
            data = data_tidy) 
summary(m2)
LMM_check(m2)
anova(m1,m2) # Significant

m3 <- glmer(spikes ~ time + N_crise_P + Type + (1|Patient),
            family = poisson(link = "log"),
            data = data_tidy) 

summary(m3)
LMM_check(m3)
anova(m2,m3) # Significant

m4 <- glmer(spikes ~ time + N_crise_P + Type + Les +(1|Patient),
            family = poisson(link = "log"),
            data = data_tidy) 

summary(m4)
LMM_check(m4)
anova(m3,m4) # Not significant

m5 <- glmer(spikes ~ time + N_crise_P + Type + localisation_temp_ext +(1|Patient),
            family = poisson(link = "log"),
            data = data_tidy) 

summary(m5)
LMM_check(m5)
anova(m3,m5) # Not significant 

m6 <- glmer(spikes ~ time + N_crise_P + Type + time:N_crise_P +(1|Patient),
            family = poisson(link = "log"),
            data = data_tidy) 

summary(m6)
LMM_check(m6)
anova(m3,m6) # Significant

m7 <- glmer(spikes ~ time + N_crise_P + Type + time:N_crise_P + time:Type + (1|Patient),
            family = poisson(link = "log"),
            data = data_tidy) 

summary(m7)
LMM_check(m7)
anova(m6,m7) # Significant


m8 <- glmer(spikes ~ time + N_crise_P + Type + time:N_crise_P + time:Type +  time*Les + (1|Patient),
            family = poisson(link = "log"),
            data = data_tidy,
            control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                   optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE))) # AIC = 10817 ; BIC = 10875.5 ; deviance = 10779
summary(m8)
LMM_check(m8)
anova(m7, m8) # Significant

m9 <- glmer(spikes ~ time + N_crise_P + Type + time:N_crise_P + time:Type + time*Les + time*localisation_temp_ext + (1|Patient),
             family = poisson(link = "log"),
             data = data_tidy,
             control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                    optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))# AIC = 10353 ; BIC = 10448 ; deviance = 10291
summary(m9)
LMM_check(m9)
anova(m8, m9) # Significant

# Summary of the models according to the anova() and the model fitting information
list_models <- list(m0 = m0, m1 = m1, m2 = m2, m3 = m3, m4 = m4, m5 = m5, 
                    m6 = m6, m7 = m7, m8 = m8, m9 = m9)
mod_fitting(list_models) # To display the models according to their fitting parameters
cross_anova_models(list_models) # To check the model selection by examining the anova of each model

m_best <- m9
summary(m_best)
LMM_check(m_best)
Anova(m_best)
r.squaredGLMM(m_best)


# post-hocs
emmeans(m_best, list(pairwise ~ time|Les), method = "tukey")
emmeans(m_best, list(pairwise ~ time| localisation_temp_ext), method = "tukey")

# Clear the environement 
rm(m0, m1, m10, m2, m3, m4, m5, m6, m7, m8, m9)
dev.off()


##------------------------------------------------------------------------
##  Influence of AEDs withdrawal on global interictal spiking activity   -
##------------------------------------------------------------------------

# data import
data_withdr <- read_excel("data_withdr.xlsx")
View(data_withdr)
data_withdr <- data_withdr%>%
  dplyr::mutate(crise = ifelse(is.na(crise), "no sz", crise)) # N_sz has many NA since not all patients experienced seizures. These NA were transformed in no_sz.

# Transformation to obtain the correct variable types
data_withdr$code <- factor(data_withdr$code)
data_withdr$onoff <- factor(data_withdr$onoff)
data_withdr$delay <- factor(data_withdr$delay)
data_withdr <- as.data.frame(data_withdr)
data_withdr$Les <- factor(data_withdr$Les)
data_withdr$Localisation <- factor(data_withdr$Localisation)
data_withdr$crise <- factor(data_withdr$crise)

# Rename variables
data_withdr <- data_withdr%>%
  rename(pat_code = code,
         N_spikes = value,
         N_day_withdr = `Delai sevrage`,
         cat_withdr_delay = delay,
         N_AED = `N molecules`,
         pres_lesion = Les,
         loc_lesion = Localisation,
         N_day_monitoring = `durée monitoring EEG`,
         sz_type = crise)


str(data_withdr)

# Data transformation to satisfy linearity --> Laplace and log-transform
data_withdr<- data_withdr%>%
  dplyr::mutate(lapl_N_spikes = N_spikes +  0.0001, # Laplacian transformation
         log_N_spikes = log(lapl_N_spikes))%>%
  dplyr::mutate(loc_temp_ext = case_when(loc_lesion == "R temporal" ~ "temporal",
                                         loc_lesion == "L temporal" ~ "temporal",
                                         TRUE ~ "extratemporal")) # Creation of a new variable groupping the localisation of the patients lesion site: either temporal or extra-temporal

data_withdr%>%
  ggplot(aes(x = log_N_spikes))+ geom_histogram()# Quite well distributed except for one bin



# Splitting of the data frame according to the onoff variable --> off1 and off2
data_withdr_off2 <- data_withdr%>%
  filter(onoff == "Off2" | onoff == "On")

data_withdr_off1 <- data_withdr%>%
  filter(onoff == "Off1"| onoff == "On")

data_withdr_off1%>%
  ggplot(aes(x=pat_code, y = log_N_spikes))+geom_point()+theme_minimal() # For each patient differences in spiking rate between on anf off

data_withdr_off2%>%
  ggplot(aes(x=pat_code, y = log_N_spikes))+geom_point()+theme_minimal() # For each patient differences in spiking rate between on anf off
# 2 patients with NA --> after removing 6h before and after the seizure --> no more data


# STATISTICAL MODELING
# Exploration of the distribution of the DV (number of spikes) and exploration of the variability of the random factors of the model
exp_plots_LMM(data_withdr, data_withdr$N_spikes, data_withdr$pat_number, data_withdr$cat_withdr_delay)  # The DV is not normally distributed (more of a log distribution since the spike variable is type double)

# Exploration of the distribution of the DV (log of number of spikes) and exploration of the variability of the random factors of the model
exp_plots_LMM(data_withdr, data_withdr$log_N_spikes, data_withdr$pat_number, data_withdr$cat_withdr_delay)  # The DV is not normally distributed (more of a log distribution since the spike variable is type double)

# The chosen model type is glmer to account for the non-linearity of the data. Since the models fits a log function, a Laplacian transformation was applied. Moreover, a log transformation of the data gave badly distributed residuals
##### Fitting of the best model for **OFF1** 
m0_glmer <- glmer(lapl_N_spikes ~ 1 + (1|pat_code),
            family = gaussian(link = "log"),
            data = data_withdr_off1) 
summary(m0_glmer)
LMM_check(m0_glmer)

m0_lmer <- lmer(log_N_spikes ~ 1 + (1|pat_code),
             REML = F,
             data = data_withdr_off1) 
summary(m0_lmer)
LMM_check(m0_lmer) # The residuals are not well distributed for a lmer model even with the log_transformation, then glmer model will be used.

m1 <- glmer(lapl_N_spikes ~ onoff + (1|pat_code),
            family = gaussian(link = "log"),
            data = data_withdr_off1,
            control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                   optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))

anova(m0_glmer, m1) # Sig
summary(m1)
LMM_check(m1)

m2 <- glmer(lapl_N_spikes ~ onoff + pres_lesion +(1|pat_code),
            family = gaussian(link = "log"),
            data = data_withdr_off1)
summary(m2)
anova(m1, m2) # Not sig
LMM_check(m2)

m3 <- glmer(lapl_N_spikes ~ onoff + loc_temp_ext +(1|pat_code),
            family = gaussian(link = "log"),
            data = data_withdr_off1)
anova(m1, m3) # Not sig
summary(m3)
LMM_check(m3)

m4 <- glmer(lapl_N_spikes ~ onoff + sz_type +(1|pat_code),
            family = gaussian(link = "log"),
            data = data_withdr_off1,
            control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                   optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
anova(m1, m4) # Not sig
summary(m4)
LMM_check(m4)

m5 <- glmer(lapl_N_spikes ~ onoff + cat_withdr_delay +(1|pat_code),
            family = gaussian(link = "log"),
            data = data_withdr_off1,
            control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                   optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
anova(m1, m5) # Not sig
summary(m5)
LMM_check(m5)

m6 <- glmer(lapl_N_spikes ~ onoff + N_AED +(1|pat_code),
            family = gaussian(link = "log"),
            data = data_withdr_off1,
            control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                   optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
anova(m1, m6) # Sig
summary(m6)
LMM_check(m6)

m7 <- glmer(lapl_N_spikes ~ onoff*pres_lesion + N_AED +(1|pat_code),
            family = gaussian(link = "log"),
            data = data_withdr_off1,
            control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                   optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
anova(m6, m7) # Sig
summary(m7)
LMM_check(m7)

m8 <- glmer(lapl_N_spikes ~ onoff*pres_lesion + N_AED + onoff*loc_temp_ext +(1|pat_code),
            family = gaussian(link = "log"),
            data = data_withdr_off1,
            control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                   optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
anova(m7, m8) # Sig
summary(m8)
LMM_check(m8)

m9 <- glmer(lapl_N_spikes ~ onoff*pres_lesion + N_AED + onoff*loc_temp_ext + onoff*sz_type + (1|pat_code),
            family = gaussian(link = "log"),
            data = data_withdr_off1,
            control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                   optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
anova(m8, m9) # Sig
summary(m9)
LMM_check(m9)

m10 <- glmer(lapl_N_spikes ~ onoff*pres_lesion + N_AED + onoff*loc_temp_ext + onoff*sz_type + onoff*cat_withdr_delay + (1|pat_code),
             family = gaussian(link = "log"),
             data = data_withdr_off1,
             control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                    optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
anova(m9, m10) # Sig
summary(m10)
LMM_check(m10)

m11 <- glmer(lapl_N_spikes ~ onoff*pres_lesion + onoff*N_AED + onoff*loc_temp_ext + onoff*sz_type + onoff*cat_withdr_delay + (1|pat_code),
             family = gaussian(link = "log"),
             data = data_withdr_off1,
             control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                    optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
anova(m10, m11) # Sig
summary(m11)
LMM_check(m11)

list_models <- list(M0 = m0_glmer, M1 = m1, M2 = m2, M3 = m3, M4 = m4, M5 = m5, M6 = m6, M7 = m7, M8 = m8, M9 = m9, M10 = m10, M11 = m11)
mod_fitting(list_models) # To display the models according to their fitting parameters
cross_anova_models(list_models) # To check the model selection by examining the anova of each model

m_best_off1 <- m11
Anova(m_best_off1)
summary(m_best_off1)
LMM_check(m_best_off1) # The postulates are acceptably respected
ICC_ranef(m_best_off1)
r.squaredGLMM(m_best_off1)


# Clear the environement 
rm(m0_lmer, m0_glmer, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11)
dev.off()

### fitting the best_model for off2

m0_glmer <- glmer(lapl_N_spikes ~ 1 + (1|pat_code),
                  family = gaussian(link = "log"),
                  data = data_withdr_off2)
summary(m0_glmer)
LMM_check(m0_glmer) # 2 Extreme residus --> to remove ? 

m1 <- glmer(lapl_N_spikes ~ onoff + (1|pat_code),
            family = gaussian(link = "log"),
            data = data_withdr_off2,
            control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                   optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))

anova(m0_glmer, m1) # Sig
summary(m1)
LMM_check(m1)

m2 <- glmer(lapl_N_spikes ~ onoff + pres_lesion +(1|pat_code),
            family = gaussian(link = "log"),
            data = data_withdr_off2,
            control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                   optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
anova(m1, m2) # Not sig
summary(m2)
LMM_check(m2)

m3 <- glmer(lapl_N_spikes ~ onoff + loc_temp_ext +(1|pat_code),
            family = gaussian(link = "log"),
            data = data_withdr_off2,
            control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                   optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
anova(m1, m3) # Not sig
summary(m3)
LMM_check(m3)

m4 <- glmer(lapl_N_spikes ~ onoff + sz_type +(1|pat_code),
            family = gaussian(link = "log"),
            data = data_withdr_off2,
            control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                   optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
anova(m1, m4) # Not sig
summary(m4)
LMM_check(m4)

m5 <- glmer(lapl_N_spikes ~ onoff + cat_withdr_delay +(1|pat_code),
            family = gaussian(link = "log"),
            data = data_withdr_off2,
            control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                   optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
anova(m1, m5) # Not sig
summary(m5)
LMM_check(m5)

m6 <- glmer(lapl_N_spikes ~ onoff + N_AED +(1|pat_code),
            family = gaussian(link = "log"),
            data = data_withdr_off2,
            control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                   optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
anova(m1, m6) # Not sig
summary(m6)
LMM_check(m6)

m7 <- glmer(lapl_N_spikes ~ onoff*pres_lesion + (1|pat_code),
            family = gaussian(link = "log"),
            data = data_withdr_off2,
            control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                   optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
anova(m1, m7) # Sig
summary(m7)
LMM_check(m7)

m8 <- glmer(lapl_N_spikes ~ onoff*pres_lesion + onoff*loc_temp_ext + (1|pat_code),
            family = gaussian(link = "log"),
            data = data_withdr_off2,
            control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                   optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
anova(m7, m8) # Sig
summary(m8)
LMM_check(m8)

m9 <- glmer(lapl_N_spikes ~ onoff*pres_lesion + onoff*loc_temp_ext + onoff*sz_type + (1|pat_code),
            family = gaussian(link = "log"), # fixed-effect model matrix is rank deficient
            data = data_withdr_off2,
            control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                   optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))

m10 <- glmer(lapl_N_spikes ~ onoff*pres_lesion + onoff*loc_temp_ext + onoff*cat_withdr_delay + (1|pat_code),
             family = gaussian(link = "log"), 
             data = data_withdr_off2,
             control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                    optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))

anova(m8, m10) # Sig
summary(m10)
LMM_check(m10)

m11 <- glmer(lapl_N_spikes ~ onoff*pres_lesion + onoff*loc_temp_ext + onoff*cat_withdr_delay + onoff*N_AED + (1|pat_code),
             family = gaussian(link = "log"), 
             data = data_withdr_off2,
             control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                    optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
anova(m10, m11) # Sig
summary(m11)
LMM_check(m11) # Some extreme residus --> to remove ?

list_models <- list(M0 = m0_glmer, M1 = m1, M2 = m2, M3 = m3, M4 = m4, M5 = m5, M6 = m6, M7 = m7, M8 = m8, M10 = m10, M11 = m11)
mod_fitting(list_models)
cross_anova_models(list_models)

m_best_off2 <- m11
LMM_check(m_best_off2)
Anova(m_best_off2)
summary(m_best_off2)
ICC_ranef(m_best_off2)
r.squaredGLMM(m_best_off2)


# Clear the environement 
rm(list_models, m0_glmer, m1, m10, m11, m2, m3, m4, m5, m6, m7, m8, m9)
dev.off()

#########################################---------END---------#########################################