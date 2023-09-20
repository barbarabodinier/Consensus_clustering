#!/bin/sh
#PBS -l walltime=20:00:00
#PBS -l select=1:ncpus=1:mem=30gb
#PBS -N real_sparcl_hclust
#PBS -J 6-9

cd /work/bbodinie/consensus_clustering/Scripts/5-public_data_weighted_clustering

module load anaconda3/personal
source activate clustering

Rscript comparison_real_weighted_sparcl.R $PBS_ARRAY_INDEX hclust

