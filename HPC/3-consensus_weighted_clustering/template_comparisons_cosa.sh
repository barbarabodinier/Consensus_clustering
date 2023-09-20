#!/bin/sh
#PBS -l walltime=100:00:00
#PBS -l select=1:ncpus=1:mem=5gb
#PBS -N simul_cosa
#PBS -q med-bio
#PBS -J 1-1010

cd /work/bbodinie/consensus_clustering/Scripts/3-consensus_weighted_clustering

module load anaconda3/personal
source activate clustering

simul_study_id={simul_study_id_input}
params_id={params_id_input}
seed=$PBS_ARRAY_INDEX
algo="hclust"

Rscript comparison_weighted_cosa.R $simul_study_id $params_id $seed $algo
