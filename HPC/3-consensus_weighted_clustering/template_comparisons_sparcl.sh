#!/bin/sh
#PBS -l walltime=30:00:00
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -N simul_sparcl
#PBS -J 1-1010
#PBS -o /dev/null
#PBS -e /dev/null

cd /work/bbodinie/consensus_clustering/Scripts/3-consensus_weighted_clustering

module load anaconda3/personal
source activate clustering

simul_study_id={simul_study_id_input}
params_id={params_id_input}
seed=$PBS_ARRAY_INDEX
#algo="hclust"
algo="kmeans"

Rscript comparison_weighted_sparcl.R $simul_study_id $params_id $seed $algo
