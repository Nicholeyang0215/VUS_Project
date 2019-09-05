# VUS_Project


This folder includes R code to conduct simulations for "A Likelihood-based Approach to Assessing Frequency of Pathogenicity among Variants of Unknown Significance (VUS) in Susceptibility Genes".


### Folders:

#### Rcode folder: all code for running simulations. 
- Fampedigree.R: script to simulate family pedigrees.

- Est_pi.R: script to estimate the proprotion of VUS that are likely pathogenic and calculate positive predictive value in simulated data. 

- Est_pitheta.R: script to esitmate the proportion of VUS that are likely pathogenic and calculate positive predictive value in simulated data with the extension of `theta`.

#### Results folder: 'Data' contains all simulated data for Figures. Figure1-4.R contain scripts to generate Figure1-4.png. 



#### Simulation Step:

1. Use `Fampredigree.R` to simulate `N` family pedigree data. 
2. Use `Est_pi.R` or `Est_pitheta.R` to estimate the proportion of VUS that are likely pathogenic in simulated family pedigree data. 
3. Repeat step1-2 2000 times for each simulation scenario. 







