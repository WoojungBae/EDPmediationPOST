#!/bin/bash
#SBATCH --job-name=EDPsim05                          # Job name
#SBATCH --array=1-500                               # At most 50 jobs at a time YMMV, make sure $
#SBATCH --cpus-per-task=1
#SBATCH --account=daniels
#SBATCH --qos=daniels-b
#SBATCH --mail-type=ALL                            # Mail events (NONE, BEGIN, END, FAIL, AL$
#SBATCH --mail-user=woojung.bae@ufl.edu    # Where to send mail
#SBATCH --mem=1G                                   # Memory limit
#SBATCH --partition=hpg-default
#SBATCH --time=2-00:00:00                          # Time limit hrs:min:sec
#SBATCH --output=outfiles/array_%A-%a_EDP05.log  # the SLURM outputs will be created i$

date; hostname; pwd

module load gcc R

#Run R script
Rscript EDP_ContinuousOutcome05.R
