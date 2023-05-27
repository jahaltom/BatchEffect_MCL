#!/bin/bash
#SBATCH -N 1
#SBATCH -p huge
#SBATCH --ntasks-per-node 40
#SBATCH -t 48:00:00
#SBATCH --mail-user=jahaltom@iastate.edu
#SBATCH --mail-type=ALL


#load required modules
conda init bash
conda  activate pyrpipe_covid
#source ~/.bashrc


# run commands
time python compute_spearman_sc.py SarsCov2_adjusted_counts.tsv   > scorrs.tsv
#time python compute_spearman_sc.py SarsCov2_Regular_counts.tsv   > scorrs.tsv

time awk -F '\t' '$3>=0.8' scorrs.tsv > scorrs_thresh_0.8.tsv
time awk -F '\t' '$3>=0.85' scorrs.tsv > scorrs_thresh_0.85.tsv
time awk -F '\t' '$3>=0.9' scorrs.tsv > scorrs_thresh_0.9.tsv

# rename annotated genes
python add_gene_names.py scorrs_thresh_0.8.tsv > scorrs_thresh_0.8_renamed.tsv
python add_gene_names.py scorrs_thresh_0.85.tsv > scorrs_thresh_0.85_renamed.tsv
python add_gene_names.py scorrs_thresh_0.9.tsv > scorrs_thresh_0.9_renamed.tsv



#load required modules
#source /ocean/projects/mcb190119p/usingh/lib/myAnacondaInstallation/etc/profile.d/conda.sh
#source ~/.bashrc
#conda activate mclpy

# run mcl
for f in $(ls *renamed.tsv)
do
  time python doMCL.py $f ${f}_mclout
done
