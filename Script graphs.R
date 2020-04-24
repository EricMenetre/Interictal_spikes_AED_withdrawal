###############################################################################################
##                                Plots relative to the paper:                               ##
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
# Eric Ménétré
# Research Assistant
# EEG and Epilepsy Unit, University Hospitals of Geneva.
# Mail:Eric.Ménétré@hcuge.ch
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
library(readxl)
library(ggpubr)
library(RColorBrewer)

##---------------------------------------------------------------
##    Changes in spiking activity during the pre-ictal period   -
##---------------------------------------------------------------

# Data import
data_N_sz <- read_excel("final data/data_N_sz.xlsx")
View(data_N_sz)

# Transformation of the data to satisfy the expectations of the tidy format
data_tidy <- data_N_sz%>%
  dplyr::select(Spikes_Baseline, Spikes_Seizure, Patient, N_crise_P, Type, Localisation, Les)%>%
  gather(key = "time", value = "spikes", -Patient, -N_crise_P, -Type, -Localisation, -Les)%>%
  mutate(spikes_w0 = spikes + 0.0001, # Laplacian transformation of the data
         log.spikes = log(spikes_w0),
         localisation_temp_ext = case_when(Localisation == "R temporal" ~ "temporal",
                                                                      Localisation == "L temporal" ~ "temporal",
                                                                      TRUE ~ "extratemporal"))

# Transformation of all the categorical variables as factor
data_tidy$Patient <- factor(data_tidy$Patient)
data_tidy$N_crise_P <- factor(data_tidy$N_crise_P)
data_tidy$Type <- factor(data_tidy$Type)
data_tidy$Localisation <- factor(data_tidy$Localisation)
data_tidy$Les <- factor(data_tidy$Les)
data_tidy$time <- factor(data_tidy$time)
data_tidy$localisation_temp_ext <- factor(data_tidy$localisation_temp_ext, levels = c("temporal", "extratemporal"))

str(data_tidy)


# Barplot of the evolution of the spike rate between baseline and pre-ictal period. The error bars represent the standard error.
p1 <- data_tidy%>%
  group_by(time)%>%
  summarise(mean_time = mean(spikes),
            sd_time = sd(spikes),
              sem_time = sd(spikes)/sqrt(nrow(data_tidy)))%>%
  mutate(newname = case_when(time == "Spikes_Baseline" ~ "Baseline",
                             time == "Spikes_Seizure" ~ "Pre-ictal"))%>%
  ggplot(aes(x = newname, y = mean_time))+
  geom_bar(stat = "identity",
           width = 0.5,
           alpha = 0.7)+
  theme_minimal()+
  geom_errorbar(aes(ymin = mean_time - sem_time,
                    ymax = mean_time + sem_time),
                width = 0.2,
                alpha = 0.8,
                size = 0.6)+
  coord_cartesian(ylim = c(50,205))+
  labs(x = "Time of the measure",
       y = "Averaged mean of spikes",
       title = "Main effect of time of measure")+
  annotate("text", x = 1.45, y = 140, label = "***", size = 9)

# Interaction time * Les
p2 <- data_tidy%>%
  group_by(time, Les)%>%
  summarise(mean_time = mean(spikes),
            sd_time = sd(spikes),
            sem_time = sd(spikes)/sqrt(nrow(data_tidy)))%>%
  mutate(newname = case_when(time == "Spikes_Baseline" ~ "Baseline",
                             time == "Spikes_Seizure" ~ "Pre-ictal"),
         Les = case_when(Les == "Les" ~ "Lesional",
                         Les == "N-Les" ~ "Non-Lesional"))%>%
  ggplot(aes(x = newname, y = mean_time, fill = Les))+
  geom_bar(stat = "identity",
           width = 0.5,
           alpha = 0.8,
           position = position_dodge())+
  theme_minimal()+
  geom_errorbar(aes(ymin = mean_time - sem_time,
                    ymax = mean_time + sem_time),
                width = 0.5,
                alpha = 0.8,
                size = 0.6,
                position = position_dodge())+
  coord_cartesian(ylim = c(50,205))+
  labs(x = "Time of the measure",
       y = "Averaged mean of spikes",
       title = "Effect of time of measure by the epilepsy type",
       fill = "Epilepsy type")+
  scale_fill_brewer(palette = "Greens")


# Interaction time * Localisation
p3 <- data_tidy%>%
  group_by(time, localisation_temp_ext)%>%
  summarise(mean_time = mean(spikes),
            sd_time = sd(spikes),
            sem_time = sd(spikes)/sqrt(nrow(data_tidy)))%>%
  mutate(newname = case_when(time == "Spikes_Baseline" ~ "Baseline",
                             time == "Spikes_Seizure" ~ "Pre-ictal"))%>%
  ggplot(aes(x = newname, y = mean_time, fill = localisation_temp_ext))+
  geom_bar(stat = "identity",
           width = 0.5,
           alpha = 0.7,
           position = position_dodge())+
  theme_minimal()+
  geom_errorbar(aes(ymin = mean_time - sem_time,
                    ymax = mean_time + sem_time),
                width = 0.5,
                alpha = 0.8,
                size = 0.6,
                position = position_dodge())+
  coord_cartesian(ylim = c(50,205))+
  labs(x = "Time of the measure",
       y = "Averaged mean of spikes",
       title = "Effect of time of measure by localisation",
       fill = "Localisation")+
  scale_fill_brewer(palette = "Purples")
  

# Interaction time*Type
p4 <- data_tidy%>%
  group_by(time, Type)%>%
  summarise(mean_time = mean(spikes),
            sd_time = sd(spikes),
              sem_time = sd(spikes)/sqrt(nrow(data_tidy)))%>%
  mutate(newname = case_when(time == "Spikes_Baseline" ~ "Baseline",
                             time == "Spikes_Seizure" ~ "Pre-ictal"))%>%
  ggplot(aes(x = newname, y = mean_time, fill = Type))+
  geom_bar(stat = "identity",
           width = 0.5,
           alpha = 0.7,
           position = position_dodge())+
  theme_minimal()+
  geom_errorbar(aes(ymin = mean_time - sem_time,
                    ymax = mean_time + sem_time),
                width = 0.5,
                alpha = 0.8,
                size = 0.6,
                position = position_dodge())+
  coord_cartesian(ylim = c(50,205))+
  labs(x = "Time of the measure",
       y = "Averaged mean of spikes",
       fill = "Seizure type",
       title = "Effect of time of measure by the seizure type")+
  scale_fill_brewer(palette = "Reds")

ggsave("First_endpoint_results.jpg", dpi = 800)
ggarrange(p1,p2,p3,p4)

# Plot of the effect of increase in spike rate depending on the seizure occurence.
data_N_sz%>%
  filter(N_crise_P <= 4)%>%
  group_by(N_crise_P)%>%
  summarize(soustr_mean = mean(soustr),
            soustr_sd = sd(soustr))%>%
  ggplot(aes(x = N_crise_P, y = soustr_mean)) + 
  geom_errorbar(aes(ymin = soustr_mean - soustr_sd, ymax = soustr_mean + soustr_sd), 
                width = 0.5, col = "seashell4")+
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = F, col = "red")+ # method = "lm",
  theme_minimal()+
  labs(x= "Seizure occurence", y = "log(spike rate pre-ictal) - log(spike rate baseline)")+
  annotate("text", x = 2.5, y = 100, label = "***", size = 9)
ggsave("Evolution seizures spike rate_means.jpg", dpi = 800)

##------------------------------------------------------------------------
##  Influence of AEDs withdrawal on global interictal spiking activity   -
##------------------------------------------------------------------------

# data import
data_withdr <- read_excel("final data/data_withdr.xlsx")
View(data_withdr)
data_withdr <- data_withdr%>%dplyr::select(-Code)%>%
  dplyr::mutate(crise = ifelse(is.na(crise), "no sz", crise)) # N_sz has many NA since not all patients experienced seizures. These NA were transformed in no_sz.

# Transformation to obtain the correct variable types
data_withdr$code <- factor(data_withdr$code)
data_withdr$Patients <- factor(data_withdr$Patients)
data_withdr$onoff <- factor(data_withdr$onoff)
data_withdr$Categories_Spikes <- factor(data_withdr$Categories_Spikes)
data_withdr$delay <- factor(data_withdr$delay)
data_withdr <- as.data.frame(data_withdr)
for(i in 11:26){data_withdr[,i] <- as.logical(data_withdr[,i])}
data_withdr$Les <- factor(data_withdr$Les)
data_withdr$Localisation <- factor(data_withdr$Localisation)
data_withdr$crise <- factor(data_withdr$crise)

# Rename variables
data_withdr <- data_withdr%>%
  rename(pat_code = code,
         pat_number = Patients,
         pct_withdr = pct.sevrage,
         N_spikes = value,
         N_spikes_cat = Categories_Spikes,
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
  dplyr::mutate(loc_temp_ext = case_when(loc_lesion == "R temporal" ~ "Temporal",
                                         loc_lesion == "L temporal" ~ "Temporal",
                                         TRUE ~ "Extratemporal")) # Creation of a new variable groupping the localisation of the patients lesion site: either temporal or extra-temporal



data_withdr%>%
  ggplot(aes(x = log_N_spikes))+ geom_histogram()# Quite well distributed except for one bin

# Splitting of the data frame according to the onoff variable --> off1 and off2
data_withdr_off2 <- data_withdr%>%
  filter(onoff == "Off2" | onoff == "On")

data_withdr_off1 <- data_withdr%>%
  filter(onoff == "Off1"| onoff == "On")


# GRAPHICAL EXPLORATION
# Effect of on vs. off1 and off2

data_withdr <- data_withdr%>%
  mutate(onoff = factor(onoff, levels = c("On", "Off1", "Off2")))%>%
  filter(log_N_spikes > -5)

p1 <- data_withdr%>%
  group_by(onoff)%>%
  summarise(MEAN = mean(lapl_N_spikes, na.rm = T),
            sem_MEAN = sd(lapl_N_spikes)/sqrt(nrow(data_withdr)))%>%
  ggplot(aes(x = onoff, y = MEAN))+
  geom_bar(width = 0.5, fill = "gray", alpha = 0.7, stat = "identity", position = position_dodge())+
  geom_errorbar(aes(ymin = MEAN - sem_MEAN, ymax = MEAN + sem_MEAN),
                width = 0.5)+
  theme_minimal()+
  labs(x = "", y = "log of mean number of spikes/h",
       title = "Effect of On vs. Off1 and Off2")

# Impact of the different epilepsy types
p2 <- data_withdr%>%
  group_by(onoff, pres_lesion)%>%
  summarise(MEAN = mean(lapl_N_spikes, na.rm = T),
            sem_MEAN = sd(lapl_N_spikes)/sqrt(nrow(data_withdr)))%>%
  mutate(pres_lesion = case_when(pres_lesion == "Les" ~ "Lesional",
                                 pres_lesion == "N-Les" ~ "Non-lesional"))%>%
  ggplot(aes(x = onoff, y = MEAN, fill = pres_lesion))+
  geom_bar(stat = "identity", position = position_dodge())+
  geom_errorbar(aes(ymin =  MEAN - sem_MEAN, ymax = MEAN + sem_MEAN), position = position_dodge())+
  theme_minimal()+
  labs(x = "", y = "Mean number of spikes/h", fill = "Epilepsy type",
       title = "Effect of the Epilepsy type")+
  scale_fill_brewer(palette = "Blues")

# Localization 
data_withdr$loc_temp_ext <- factor(data_withdr$loc_temp_ext, levels = c("Temporal", "Extratemporal"))
p3 <- data_withdr%>%
  group_by(onoff, loc_temp_ext)%>%
  summarise(MEAN = mean(lapl_N_spikes, na.rm = T),
            sem_MEAN = sd(lapl_N_spikes)/sqrt(nrow(data_withdr)))%>%
  ggplot(aes(x = onoff, y = MEAN, fill = loc_temp_ext))+
  geom_bar(stat = "identity", position = position_dodge())+
  geom_errorbar(aes(ymin =  MEAN - sem_MEAN, ymax = MEAN + sem_MEAN), position = position_dodge())+
  theme_minimal()+
  labs(x = "", y = "Mean number of spikes/h", fill = "Localisation",
       title = "Effect of localisation of epileptic foyer")+
  scale_fill_brewer(palette = "Reds")

# Type of seizure
p4 <- data_withdr%>%
  mutate(sz_type = as.character(sz_type),
         sz_type = ifelse(sz_type == "no sz", "No seizure", sz_type),
         sz_type = factor(sz_type))%>%
  group_by(onoff, sz_type)%>%
  summarise(MEAN = mean(lapl_N_spikes, na.rm = T),
            sem_MEAN = sd(lapl_N_spikes)/sqrt(nrow(data_withdr)))%>%
  ggplot(aes(x = onoff, y = MEAN, fill = sz_type))+
  geom_bar(stat = "identity", position = position_dodge())+
  geom_errorbar(aes(ymin =  MEAN - sem_MEAN, ymax = MEAN + sem_MEAN), position = position_dodge())+
  theme_minimal()+
  labs(x = "", y = "Mean number of spikes/h", fill = "Types of seizure",
       title = "Effect of the type of seizure")+
  scale_fill_brewer(palette = "Greens")

# Category of withdrawal delay
p5 <- data_withdr%>%
  dplyr::mutate(cat_withdr_delay = case_when(cat_withdr_delay == 1 ~ "2 to 4 days",
                                             cat_withdr_delay == 2 ~ "5 to 8 days",
                                             cat_withdr_delay == 3 ~ "9 to 11 days"),
                cat_withdr_delay = factor(cat_withdr_delay, levels = c("2 to 4 days", "5 to 8 days", "9 to 11 days")))%>%
  group_by(onoff, cat_withdr_delay)%>%
  summarise(MEAN = mean(lapl_N_spikes, na.rm = T),
            sem_MEAN = sd(lapl_N_spikes)/sqrt(nrow(data_withdr)))%>%
  ggplot(aes(x = onoff, y = MEAN, fill = cat_withdr_delay))+
  geom_bar(stat = "identity", position = position_dodge())+
  geom_errorbar(aes(ymin =  MEAN - sem_MEAN, ymax = MEAN + sem_MEAN), position = position_dodge())+
  theme_minimal()+
  labs(x = "", y = "Mean number of spikes/h", fill = "Withdrawal delays",
       title = "Effect of the Withdrawal delay")+
  scale_fill_brewer(palette = "YlGnBu")

# Number of AED
p6 <- data_withdr%>%
  dplyr::mutate(N_AED = factor(N_AED))%>%
  group_by(onoff, N_AED)%>%
  summarise(MEAN = mean(lapl_N_spikes, na.rm = T),
            sem_MEAN = sd(lapl_N_spikes)/sqrt(nrow(data_withdr)))%>%
  ggplot(aes(x = onoff, y = MEAN, fill = N_AED))+
  geom_bar(stat = "identity", position = position_dodge())+
  geom_errorbar(aes(ymin =  MEAN - sem_MEAN, ymax = MEAN + sem_MEAN), position = position_dodge())+
  theme_minimal()+
  labs(x = "", y = "log of mean number of spikes/h", fill = "Number of AED",
       title = "Effect of the number of AED")+
  scale_fill_brewer(palette = "Purples")

ggarrange(p1, p2, p3, p4, p5, p6)
ggsave("Second endpoint results.jpg", dpi = 800)
rm(p1, p2, p3, p4, p5, p6)

# DESCRIPTIVE STATISTICS

# Exploration des valeurs

data_withdr%>%
  group_by(onoff)%>%
  summarise(mean_raw = round(mean(N_spikes, na.rm = T),2), 
            sd_raw = round(sd(N_spikes, na.rm = T),2),
            mean_log = round(mean(log_N_spikes, na.rm = T),2),
            sd_log = round(sd(log_N_spikes, na.rm = T),2))


data_withdr%>%
  mutate(pres_lesion = case_when(pres_lesion == "Les" ~ "Lesional",
                                 pres_lesion == "N-Les" ~ "Non-lesional"))%>%
  group_by(onoff, pres_lesion)%>%
  summarise(mean_raw = round(mean(N_spikes, na.rm = T),2), 
            sd_raw = round(sd(N_spikes, na.rm = T),2),
            mean_log = round(mean(log_N_spikes, na.rm = T),2),
            sd_log = round(sd(log_N_spikes, na.rm = T),2))

data_withdr%>%
  group_by(onoff, loc_temp_ext)%>%
  summarise(mean_raw = round(mean(N_spikes, na.rm = T),2), 
            sd_raw = round(sd(N_spikes, na.rm = T),2),
            mean_log = round(mean(log_N_spikes, na.rm = T),2),
            sd_log = round(sd(log_N_spikes, na.rm = T),2))


data_withdr%>%
  filter(!is.na(sz_type))%>%
  group_by(onoff, sz_type)%>%
  summarise(mean_raw = round(mean(N_spikes, na.rm = T),2), 
            sd_raw = round(sd(N_spikes, na.rm = T),2),
            mean_log = round(mean(log_N_spikes, na.rm = T),2),
            sd_log = round(sd(log_N_spikes, na.rm = T),2))

