# VUS_Project


This folder includes R code to conduct simulations for "A Likelihood-based Approach to Assessing Frequency of Pathogenicity among Variants of Unknown Significance (VUS) in Susceptibility Genes".

- R code: 
1. Fampedigree.R: script to simulate family pedigrees.
2. Est_pi.R: script to estimate the proprotion of VUS that are likely pathogenic in simulated data. 
3. Est_pitheta.R: script to esitmate the proportion of VUS that are likely pathogenic in simulated data with the extension of `theta`.


- Simulation Step:

1. Use `Fampredigree.R` to simulate `N` family pedigree data. 
2. Use `Est_pi.R` or `Est_pitheta.R` to estimate the proportion of VUS that are likely pathogenic in simulated family pedigree data. 
3. Repeat step1-2 2000 times for each simulation scenario. 







