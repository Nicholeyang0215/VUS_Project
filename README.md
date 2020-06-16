# VUS_Project


This folder includes R code to conduct simulations for "A Likelihood-based Approach to Assessing Frequency of Pathogenicity among Variants of Unknown Significance (VUS) in Susceptibility Genes".


### Folders:

#### 1. Rcode folder: all code for running simulations. 
- sim_fam_varied.R: script to simulate family pedigrees based on different family structures in the main text.
- sim_fam_constant.R: script to simulate family pedigrees based on a constant family structure.


- estimate_pi.R: script to estimate the proprotion of VUS that are likely pathogenic and calculate positive predictive value in simulated data. 

- estimate_pi_with_theta.R: script to esitmate the proportion of VUS that are likely pathogenic and calculate positive predictive value in simulated data with the extension of `theta`.

#### 2. Results folder: 'Data' contains simulated data for Figures 1-6. It also contains scripts for generating Figure1-6. 




### Simulation Step:

1. Use `sim_fam_varied.R` to simulate family pedigree data of size `N`. Under each sample size, replicates are named with `script_num.RData`.  
2. Use `estimate_pi.R` or `estimate_pi_with_theta.R` to estimate the proportion of VUS that are likely pathogenic and the positive predictive value in simulated family pedigree data. 
3. Repeat step1-2 for all simulated cohort replicates. 








