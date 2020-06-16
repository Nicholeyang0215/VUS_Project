
# Description: This folder contains estimation results based on simulated cohort data. 

1. Those .csv files are summarized results for generating Figure 2-3 in the main text. The setting of sample size = 1 million
contain 50 replicates, while all other scenarios contain 2000 estimation replicates. 

2. The `Raw_estimation_results` folder contains 3 examples of estimation results for Figure 2(a). Those are replicates with
script_num=1 under the simulated cohorts of '1k_fam1.RData', '5k_fam1.RData' and '1w_fam1.RData'.

3. The sbatch command for get all the raw estimation results is as below: 

sbatch --array=1-2500 estimate_pi.sh
