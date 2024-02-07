#!/bin/bash 

#SBATCH --job-name=SSIZE06K
#SBATCH --partition=short
#SBATCH --nodes=1 
#SBATCH --mem=250G 
#SBATCH --time=04:00:00
#SBATCH --output=%j.output 
#SBATCH --error=%j.error

date
module load R/3.6.2 
R CMD BATCH cal.metrics.sim_circuits.R
date

