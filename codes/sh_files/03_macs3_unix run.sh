#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p pq #highmem # partition to use
#SBATCH --time=10:00:00 # maximum walltime
#SBATCH -A xxxx # research project to submit under
#SBATCH --nodes=1 # specify number of nodesq p
#SBATCH --exclusive
#SBATCH --job-name=macs
#SBATCH --mem=100G #400G
#SBATCH --error=macs_%A_%a.err

mamba init
mamba activate /lustre/home/mc900/RP_T121869/.conda/envs/myenvR

Rscript ../03_2_macs3_unix run.R
