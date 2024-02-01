# Script: General pipeline to calculate classical summary indicators of cortisol dynamics in the STRESS-EU database (www.stressdatabase.eu)
# Authors: Milou Sep, Laura de Nooij and Jonathan Posthuma
# Contact: m.s.c.sep(at)amsterdamumc.nl

# Input: central stored anonymous data (note, see www.stressdatabase.eu for the data access procedure)
# participants.csv contains all cortisol values and time points of the studies in the STRESS-EU database
# studies.csv contains all meta-information on the included studies in the STRESS-EU database

# Environment preparation -------------------------------------------------
rm(list=ls())

# Load packages
library(readr)
library(dplyr)
library(tidyr)
library(purrr)

# Load data
participants <- read_csv("data/participants.csv", na =c("NA", "not_available", "not_measured", "not_applicable", "not available",""))
studies <- read_csv("data/studies.csv",  na =c("NA", "not_available", "not_measured", "not_applicable", "not available",""))

# Data preparation --------------------------------------------------------
# Recode missing values
participants %>% 
  mutate_all(~gsub('not_available|not_applicable|NA|not_Available|not_availble|not_measured', NA, .)) -> participants_recoded
# Remove rows without cortisol values
participants_recoded %>%
  filter(!is.na(cortisol_timepoint) &  !is.na(cortisol_value)) -> participants_cortisol
# Check uniformity cortisol units
unique(participants_cortisol$cortisol_unit) # Note, correct to nmol/l if needed

# Baseline Correction -----------------------------------------------------
baseline_correction <- function(t,v){
  # t = timepoints of all data samples, calibrated to t=0 as onset stressor (vector)
  # v = values of cortisol (nmol/L) corresponding to timepoints (vector)
  
  # Change character vectors to numeric vector
  strsplit(t, ",") %>% unlist(.) %>% gsub("NA",NA,.) %>% as.numeric(.)->t
  strsplit(v, ",") %>% unlist(.) %>% gsub("NA",NA,.) %>% as.numeric(.)->v
  
  # Identify location baseline values (i.e. between 20min before stressor onset up to 2 min after stressor onset)
  
  # Drop values more then 20min prior stressor onset if applicable
  if(any(t < -20)){
    which(t< -20) -> index_redundant_baseline_time_points
    v[-c(index_redundant_baseline_time_points)]->v
    t[-c(index_redundant_baseline_time_points)]->t
  }
  
  # Select baseline values to use (between -20min prior and 2 min after stressor onset) to use (Note, values <-20 or > 2 are not used for the baseline correction) 
  which(t %in% -20:2) -> index_baseline_time_points
  
  if(length(index_baseline_time_points) == 0) { # Correction for studies without baseline measures, make missing explicit (e.g. add t0 = NA)
    print('no baseline values')
    c(NA, v)->v
    c(0, t)->t
    
  } else if ( # Correction for 1 baseline value below 0: if one baseline value t<0, then also set to t=0 as baseline
    length(index_baseline_time_points) == 1){
    if(t[index_baseline_time_points] < 0){
      t[index_baseline_time_points] <- 0
    } 
    
  } else if (length(index_baseline_time_points) > 1 ) { # Correction for multiple baseline measures (Take average of multiple baseline values and assign to time point 0)
    v[index_baseline_time_points]-> baseline_values
    c(mean(baseline_values), v[-index_baseline_time_points])->v
    c(0, t[-index_baseline_time_points])->t
  } 
  
  # Collect variables: n_baseline, t and v
  out<-list(n_t=length(index_baseline_time_points),t=base::paste(t, collapse = ","),v=base::paste(v, collapse = ","))
  
  # Return output
  return(out)
}

# Call function baseline correction
participants_cortisol %>%
  rowwise() %>%
  mutate(
    base_corr=  map2(cortisol_timepoint, cortisol_value, baseline_correction),
    .keep='all'
  ) %>% 
  mutate(
    n_baseline = unlist(base_corr)[1],
    cortisol_timepoint = unlist(base_corr)[2],
    cortisol_value = unlist(base_corr)[3]
  ) %>% select(-base_corr)->participants_baseline_corrected

# Summary cortisol indicators ---------------------------------------------

# Information on formula's for the summary cortisol indicators 

# [Pruessner, 2003](https://www.sciencedirect.com/science/article/pii/S0306453002001087#FD6) provide two formula's for the calculation of the ‘Area under the curve with respect to increase’ (AUCi) and ‘Area under the curve with respect to ground’ (AUCg) and reflect on the differences in meaning between these formula's.
# - The information needed in order to calculate the formula consists of the values themselves (v) and the time distance between the values (t)
# - [Khoury, 2015](https://www.sciencedirect.com/science/article/pii/S2352289515000272) confirm that the AUCi and AUCg represent different aspects of the stress response (not based on cortisol!). They also provide information on the interpretation AUCg and AUCi and the rationale behind the selection of these two sumscores.
# - The formula from Pruessner, et al (2003) for the __AUCg__ is $AUCg = ( ( (v2 + v1)*t1 )/2 ) + ( ( (v3 + v2)*t2 )/2 ) + ( ( (v4 + v3)*t3 )/2 ) + ( ( (v5 + v4)*t4 )/2 ) + ( ( (v6 + v5)*t5 )/2 ) + etc.$ In which v1 to v6 denote the single values, and t1 to t5 denote the time interval between the values.
# - This formula is transformed into $AUCg = SUM( ( (v[i+1] + v[i])*t[i] )/2 )$ and used in a for-loop.
# - The formula for the __AUCi__ can be derived from the formula for AUCG, since it is identical to AUCG except for the removal of the area between ground and the first measure (baseline) for all time points (Pruessner, 2003): $AUCi = AUCg - (v[1] * SUM(time.intervals))$
# function adapted from https://github.com/mscsep/SAM_explore/blob/main/R/SAMexplore_AUC_calculations.Rmd

# other summary measures: Khoury JE, Gonzalez A, Levitan RD, Pruessner JC, Chopra K, Basile VS, Masellis M, Goodwill A, Atkinson L. Summary cortisol reactivity indicators: Interrelations and meaning. Neurobiol Stress. 2015 Apr 30;2:34-43. doi: 10.1016/j.ynstr.2015.04.002. PMID: 26844238; PMCID: PMC4721456.
# - Peak reactivity (khoury, 2015: Change in cortisol between the baseline and peak values.= peak (max value after baseline) - baseline)
# - Slope of the line between baseline and peak cortisol value. slope_reactivity = (peak - baseline)/ (time_peak - time_baseline).
# - Slope of line between peak and last cortisol value (max 2h after stressor): slope_recovery = (peak - last_value)/(time_last_value - time_peak)

summary_indicators <- function(t, v){
  # t = timepoints of all data samples, calibrated to t=0 as onset stressor (vector)
  # v = values of cortisol (nmol/L) corresponding to timepoints (vector)
  # Note, this function requires 1 baseline value; use function 'baseline_correction' in this script first
  
  # 1A) Change character vectors to numeric vector
  strsplit(t, ",") %>% unlist(.) %>% gsub("NA",NA,.) %>% as.numeric(.)->t
  strsplit(v, ",") %>% unlist(.) %>% gsub("NA",NA,.) %>% as.numeric(.)->v
  
  # 1B) Check for if there are duplicate timepoints, if present remove latest from t and v
  if (any(duplicated(t))){
    print("duplicate timepoints")
    which(t == t[anyDuplicated(t)]) -> index_duplicated_timepoints
    t = t[-max(index_duplicated_timepoints, na.rm = T)]
    v = v[-max(index_duplicated_timepoints, na.rm = T)]
  }
  
  # 2A) Only proceed with studies that have a baseline value (e.g. t=0 is not NA) and more than 1 time point; AND if vector t and v have the same length
  if ((length(v) == length(t)) & 
      (!is.na(v[which(t == 0)]) & length(t)>1)){
    
    # Set limits for calculations (up to 2h post stressor onset)
    which(t < 120) -> index_max_2h_values
    t[index_max_2h_values] -> t_max2h
    v[index_max_2h_values] -> v_max2h
    
    # 3) Peak reactivity (based on Khoury, 2015): Change in cortisol between the baseline and peak values.
    # 3A) find index of baseline values (NB, should be 1 after baseline correction)
    which(t_max2h<=0) -> index_baseline_time_points
    # 3B) find index of peak value within 15-50 minutes post baseline
    which(t_max2h %in% c(15:50)) -> index_range_peak_value
    # Set peak and slope reactivity to NA if no values within peak range (15-50 min.)
    if (length(index_range_peak_value) == 0){
      peak_reactivity = NA
      slope_reactivity = NA
      slope_recovery = NA
      # if values within peak range (15-50 min)  
    } else if (length(index_range_peak_value) !=0){
      which(v==max(v_max2h[index_range_peak_value], na.rm = T)) -> index_peak_value # NB refers to indices of v
      # 3B.1) Correction for multiple identical peak values: select the first value (~ start peak)
      if (length(index_peak_value) >1){index_peak_value = min(index_peak_value, na.rm = T)}
      # 3C) peak_reactivity = peak (max value after baseline) - baseline). Calculate peak reactivity if baseline value is not identical to peak value, else NA
      peak_reactivity = v_max2h[index_peak_value] - v_max2h[index_baseline_time_points]
      
      # 4) Slope reactivity (based on Khoury, 2015): Slope of the line between baseline and peak cortisol value. slope_reactivity = (peak - baseline)/ (time_peak - time_baseline)
      slope_reactivity = (v_max2h[index_peak_value] - v_max2h[index_baseline_time_points]) / (t_max2h[index_peak_value] - t_max2h[index_baseline_time_points])
      
      # 5) Slope recovery (following Khoury, 2015 formula for slope reactivity): Slope of the line between peak and last cortisol value
      # 5A) find last time point
      which(t_max2h==max(t_max2h, na.rm = T)) -> index_last_value # waaorm max?
      # 5B) slope_recovery = (peak - last value) / (time_last value - time_peak )
      slope_recovery = ( v_max2h[index_last_value]- v_max2h[index_peak_value]) / (t_max2h[index_last_value] - t_max2h[index_peak_value])
    }
    
    # 6) AUCg & AUCi calculations (based on formula in Pruessner, 2003)
    # 6A) calculate time intervals
    t_int=diff(t_max2h) # diff is used to calculated intervals between the time points
    # 6B) Calculate AUC per interval
    list_AUC<-list()
    for(i in 1:(length(t_int))){
      (((v_max2h[i+1] + v_max2h[i]) * t_int[i] ) / 2) -> list_AUC[[i]]
      }
    # 6C) change list to vector
    unlist(list_AUC)->vector_AUC
    # 6D) AUCg: sum AUCs of each time interval
    sum(vector_AUC, na.rm = F)->AUCg
    # AUCi: remove area between ground and the first measurement (v1), for all time points (=sum(time_intervals) from AUCg (Pruessner, 2003).
    # 6E) calculate baseline:
    if(length(index_baseline_time_points) == 1){
      baseline = v_max2h[index_baseline_time_points] * sum(t_int) 
      # 6F) subtract baseline from AUCg
      AUCi<- AUCg - baseline 
    } else if (length(index_baseline_time_points) > 1){
      AUCi = NA
    }
    
    # 7) Count number of time points used for summary indicator calculations
    length(t_max2h)-> n_timepoints
    
    # 2B) Set summary indicator values to missing, for studies with no baseline and less than 1 time point
  } else { 
    AUCg = NA
    AUCi = NA
    peak_reactivity = NA
    slope_reactivity = NA
    slope_recovery = NA
    n_timepoints = NA
  }
  
  # 8) Collect summary indicator variables
  out<-list(AUCg, AUCi, peak_reactivity, slope_reactivity, slope_recovery, n_timepoints)
  return(out)
}          

# Call function summary indicators
participants_baseline_corrected %>%
  rowwise() %>%  # compute summary indicators per row
  mutate(
    AUCg = summary_indicators(cortisol_timepoint, cortisol_value)[[1]],
    AUCi = summary_indicators(cortisol_timepoint, cortisol_value)[[2]],
    peak_reactivity = summary_indicators(cortisol_timepoint, cortisol_value)[[3]],
    slope_reactivity = summary_indicators(cortisol_timepoint, cortisol_value)[[4]],
    slope_recovery = summary_indicators(cortisol_timepoint, cortisol_value)[[5]],
    n_t_summary_indicators= summary_indicators(cortisol_timepoint, cortisol_value)[[6]]
  ) -> participants_summary_cortisol

# Add study information ---------------------------------------------------

# Get study_id variable names
names(participants_summary_cortisol) #"studyID"
names(studies) #"_id"

# Check for missing study_ids (should be non)
participants_summary_cortisol %>% filter(is.na(studyId))
studies %>% filter(is.na("_id"))

# rename study_id variables
participants_summary_cortisol = rename(participants_summary_cortisol, study_id = studyId)
studies = rename(studies, study_id = "_id")

# Join participants and study information
full_join(participants_summary_cortisol, studies, by= 'study_id')->df_merged
summary(df_merged)

# Format joined dataset
df_merged %>%
  mutate_at(c("gender", "contraceptive", "acute_stress_test", "stress_control_condition", 
              "diagnosis", grep("_measured", names(.), value=T), "biobank_stored", 
              "questionnaires_vas_yes_no", "time_of_day" ), as.factor) %>%
  mutate_at(c("age_years", "bmi"), as.numeric) %>% # Note, missings introduced because in some rows bmi has the character value: 'normal'
  droplevels()->df_formated
summary(df_formated)

# Save joined dataset -----------------------------------------------------
saveRDS(df_formated, 'processed_data/df_sum_cort.rds')