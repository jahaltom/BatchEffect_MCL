#!/bin/bash
#SBATCH -N 1
#SBATCH -p defq
#SBATCH --ntasks-per-node 8
#SBATCH -t 2:00:00
#SBATCH --mail-user=haltomj@chop.edu
#SBATCH --mail-type=ALL


module load R

Rscript GoEnrich.r
