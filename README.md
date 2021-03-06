# EpistasisStats - Data and scripts for analyzing Epistasis

### Clone of https://github.com/YtrebergPatelLab/EpistasisStats, I'm the original author. 


## General overview

### Analysis

This folder contains the scripts used for analysis. This includes:

1. Loglikelihood ratio analysis: Generates simulated set of lambas by greatest likelihood
2. Model building/selection: Performs model selection using data curated from the SKEMPIv2.0 and ProTherm4 datasets
3. Validation: Performs a leave-10-percent-systems-out validation process on the aforementioned model selection procedure. 

### Data

Contains the data used for input and generated by the analysis scripts

1. loglikedata: Data generated by the loglikelihood ratio test, used to generate figure 3
2. Mappings: Amino acid data to map the single letter code to attributes, like charge.
3. Model_build_results: Resulting best model and delta r-squared data. 
4. processed: The processed datasets from SKEMPIv2.0 and ProTherm4
5. raw: The raw datasets from SKEMPIv2.0 and ProTherm4
6. validation_results: The results from the validation procedure
7. validation_results_parsed: The split and parsed validation results

### Figures

Contains figures present in the manuscript and the scripts used to generate them. Self explanatory, for more details see the included README.md
