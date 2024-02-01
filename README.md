# General pipeline to calculate classical summary indicators of cortisol dynamics in the STRESS-EU database (www.stressdatabase.eu)
Authors: Milou Sep, Laura de Nooij and Jonathan Posthuma
contact: m.s.c.sep(at)amsterdamumc.nl

# Index

_ `README.md`: an overview of the project
|___ `data`: folder for raw data files (note, see www.stressdatabase.eu for the data access procedure)
|___ `r`: folder with r-scripts
|___ `processed_data`: folder for processed data files

R-script `stresseu_summary_indicators_cortisol_dynamics.r` contains:
* Baseline correction:
  - for studies with multiple baseline values: Replace baseline values by average
  - for studies without baseline values: Makes NA explicit on T0

* Calculation of classical summary indicators of cortisol dynamics:
  - Peak reactivity: peak - baseline
  - Slope reactivity: slope baseline to peak
  - Slope recovery: slope peak to last value
  - Area under the curve with respect to ground (AUCg)
  - Area under the curve with respect to increase (AUCi) 

For background and formula's for the summary cortisol indicators see [Khoury, 2015](https://www.sciencedirect.com/science/article/pii/S2352289515000272) and [Pruessner, 2003](https://www.sciencedirect.com/science/article/pii/S0306453002001087#FD6).