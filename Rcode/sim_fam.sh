#!/bin/bash

#SBATCH -n 1 # cores
#SBATCH -N 1 # nodes
#SBATCH -p serial_requeue # Partition
#SBATCH -t 0-00:30 # (D-HH:MM)
#SBATCH --mem 6500 # Memory request
#SBATCH -o brcapro.out # Standard output
#SBATCH -e brcapro.err # Standard error
#SBATCH --mail-type=END

module load R/3.5.1-fasrc01

seed=${SLURM_ARRAY_TASK_ID}
argString="--args "$seed

R CMD BATCH --quiet --no-restore --no-save "$argString" sim_fam_varied.R
