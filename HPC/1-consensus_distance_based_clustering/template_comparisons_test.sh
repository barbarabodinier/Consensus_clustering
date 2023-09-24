#!/bin/sh
#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -N simul_comp
#PBS -q med-bio

cd /work/bbodinie/consensus_clustering/Scripts/1-consensus_distance_based_clustering

module load anaconda3/personal
source activate clustering

simul_study_id=1
params_id=1
seed=1
algo="pam"

Rscript comparison.R $simul_study_id $params_id $seed $algo
