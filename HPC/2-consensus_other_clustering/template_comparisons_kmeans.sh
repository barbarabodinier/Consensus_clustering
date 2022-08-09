#!/bin/sh
#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=1:mem=5gb
#PBS -N simul_comp_kmeans
#PBS -q med-bio
#PBS -J 1-1000

cd /work/bbodinie/consensus_clustering/Scripts/2-consensus_other_clustering
module load anaconda3/personal

simul_study_id={simul_study_id_input}
params_id={params_id_input}
seed=$PBS_ARRAY_INDEX

Rscript comparison_kmeans.R $simul_study_id $params_id $seed
