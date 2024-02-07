#!/bin/bash 

#SBATCH --partition=short
#SBATCH --job-name=simPara
#SBATCH --time=24:00:00
#SBATCH --nodes=2 
#SBATCH --ntasks-per-node=28 
#SBATCH --mem=250G 
#SBATCH --output=%j.output 
#SBATCH --error=%j.error

date
module load R/3.6.2 
R CMD BATCH sim.circuits.para.savedbyfile.R
date

