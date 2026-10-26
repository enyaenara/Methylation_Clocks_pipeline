# Scripts
Scripts to accompany the manuscript: 
"Fanconi anemia is associated with signatures of accelerated biological aging".

This directory contains R scripts for data processing, epigenetic clock calculation, statistical analysis, and plots.

## Pipeline Overview:

- **01_Process_FAvsControls.R**: 
  Performs quality control on IDAT files, sample/probe filtering, normalization, and data harmonization across EPIC v1 and v2 datasets for the main analysis.

- **02_Clocks_FAvsControls.R** 
  Calculates epigenetic age and its acceleration across multiple clocks and generates regression and boxplot comparisons.

- **03_Z-score_calculation.R** 
  Calculates Z-scores for EAA relative to the control group and adjust biological variables

- **04_Individual_analysis.R** 
  Generates plots comparing chronological age and EAA across tissues, alongside age-group comparisons.
  
 - **05_Paired_samples_analysis.R** 
  Performs paired statistical analyses comparing blood and bone marrow samples from the same FA patient.
  
- **06_Process_FA_WS_HGPSvsControls.R**: 
  Performs quality control on IDAT files, sample/probe filtering, normalization, and data harmonization across EPIC v1 and v2 datasets for the supplementary analysis.

- **07_Episignature_FA_WS_HGPSvsControls.R** 
 Analyzes the Fanconi anemia DNA methylation signature (MDS, heatmaps) and identifies Differentially Methylated Regions (DMRs) using DMRcate.

- **08_Analysis_FA_WS_HGPSvsControls.R** 
  Performs comparative analyses incorporating progeroid syndrome individuals (WS, HGPS, NCPL) alongside Fanconi anemia and controls.