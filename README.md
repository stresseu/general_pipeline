# General pipeline to calculate classical summary indicators of cortisol dynamics in the STRESS-EU database (www.stressdatabase.eu)
- Authors: Milou Sep, Laura de Nooij and Jonathan Posthuma
- contact: m.s.c.sep(at)amsterdamumc.nl

# Index
- `README.md`: an overview of the project
- `data`: folder for raw data files (note, see www.stressdatabase.eu for the data access procedure)
- `r`: folder with r-scripts
- `processed_data`: folder for processed data files

## R-script `stresseu_summary_indicators_cortisol_dynamics.r` contains:
* Input: `data/participants.csv` and `data/studies.csv` (see www.stressdatabase.eu for the data access procedure)
* Output: `processed_data/df_sum_cort.rds`

* Baseline correction:
  - Timepoints: values between 20 minutes before and 2 minutes after stressor onset are considered baseline
  - for studies with multiple baseline values: replace baseline values by average
  - for studies without baseline values: makes NA explicit on T0

* Calculation of classical summary indicators of cortisol dynamics:
  - Timepoints: values more then 2h after stressor onset are excluded, peak values are determined between 15 and 50 minutes after stressor onset
  - Peak reactivity: peak - baseline
  - Slope reactivity: slope baseline to peak
  - Slope recovery: slope peak to last value
  - Area under the curve with respect to ground (AUCg)
  - Area under the curve with respect to increase (AUCi) 

# References
For background and formula's for the summary cortisol indicators see [Khoury, 2015](https://www.sciencedirect.com/science/article/pii/S2352289515000272) and [Pruessner, 2003](https://www.sciencedirect.com/science/article/pii/S0306453002001087#FD6).
