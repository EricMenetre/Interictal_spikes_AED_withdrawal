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
library(glmmTMB)

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
  mutate(localisation_temp_ext = case_when(Localisation == "R temporal" ~ "temporal",
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

# STATISTICAL MODELING

mod_endpoint_1 <- glmmTMB(spikes ~ time:N_crise_P + time:Type + time*localisation_temp_ext + time*Les + (1|Patient),
              data = data_tidy,
              ziformula = ~ time:N_crise_P + time:Type + time*localisation_temp_ext + time*Les,
              family = truncated_poisson)
summary(mod_endpoint_1)
Anova(mod_endpoint_1)


# post-hocs
emmeans(mod_endpoint_1, list(pairwise ~ time| localisation_temp_ext), method = "tukey")


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
# The chosen model type is glmer to account for the non-linearity of the data. Since the models fits a log function, a Laplacian transformation was applied. Moreover, a log transformation of the data gave badly distributed residuals
##### OFF 1
m_off_1 <- glmer(lapl_N_spikes ~ onoff*pres_lesion + onoff*N_AED + onoff*loc_temp_ext + onoff*sz_type + onoff*cat_withdr_delay + (1|pat_code),
             family = gaussian(link = "log"),
             data = data_withdr_off1,
             control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                    optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE))) 


Anova(m_off_1)
r.squaredGLMM(m_off_1)

### fitting the best_model for off2

m_off_2 <- glmer(lapl_N_spikes ~ onoff*pres_lesion + onoff*N_AED + onoff*loc_temp_ext  + onoff*cat_withdr_delay +  (1|pat_code),
                 family = gaussian(link = "log"),
                 data = data_withdr_off2,
                 control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                        optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE))) 
# sz_type could not be included in that model since the lower number of data points

Anova(m_off_2)
r.squaredGLMM(m_off_2)

#########################################---------END---------#########################################