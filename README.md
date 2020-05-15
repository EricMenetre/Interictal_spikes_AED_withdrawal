# Interictal_spikes_AED_withdrawal

## Data dictionary

### data_tidy (transformation of data_N_sz)
* Patient: code of each patient, coded as factor
* N_crise_P: number of seizure occurence per patient, coded as numeric
* Type: type of seizure, either GTCS, Focal or No seizure, coded as factor
* Localisation: localisation of the epileptic foyer by lobe (temporal - extratemporal), coded as factor
* Les: type of epilepsy, either lesional of non-lesional
* time: time of measure of the spiking activity, either duting baseline or in the pre-ictal period, coded as factor
* spikes: number of epileptic spikes detected, coded as numerical
* spikes_w0: same as previous one but with a Laplacian transformation (add of 0.0001 to each value to be able to apply log transformation), coded as numerical
* log.spikes: log values of the spikes_w0, coded as numerical
* location_temp_ext: clustering of the localisation information into temporal and extratemporal regions, coded as factor

### data_withdr / data_withdr_off1 / data_withdr_off2 (after transformations)
* pat_code: code of each patient, coded as factor
* onoff: time of measure of the spiking activity. Either during ON (full dose medication) vs. OFF1 (first day of maximum withdrawal except 3h before and after the seizure if a seizure occurs) and OFF2 (first day of maximum withdrawal except 6h before and after the seizure if a seizure occurs). This distinction was performed to get the results comparable with Goncharova and colleagues. This variable was coded as factor
* % AEDs: percentage of the full dose of anti-epileptic drug, coded as numerical
* N_spikes: average number of spikes per hour during the 24h either in ON or OFF, coded as numeric
* N_day_withdr: number of days until maximum withdrawal, coded as numeric
* cat_N_day_withdr: categories based on the number of days to reach the maximum withdrawal (fast, medium or slow), coded as factor
* N_AED: number of anti-epileptic drugs taken by a patient, coded as numeric
* pres_lesion: epilepsy type, either lesional or non-lesional, coded as factor
* loc_lesion: localisation of the epileptic foyer, coded as factor
* N_day_monitoring: number of days of monitoring, coded as numeric 
* sz_type: type of seizure, either focal, tonic-clonic bilateral or no seizure, coded as factor
* lapl_N_spikes: number of spikes with a Laplacian transformation (add of 0.0001 to each value), coded as numeric
* log_N_spikes: log of the laplacian number of spikes, coded as numeric
* loc_temp_ext: localisation, either temporal or extratemporal, coded as factor

## Contact information

Eric.Menetre@unige.ch

Eric.Menetre@hcuge.ch

