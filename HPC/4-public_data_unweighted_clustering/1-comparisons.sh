#!/bin/sh
#PBS -l walltime=50:00:00
#PBS -l select=1:ncpus=1:mem=20gb
#PBS -N real_unweighted
#PBS -J 1-5

cd /work/bbodinie/consensus_clustering/Scripts/4-public_data_unweighted_clustering

module load anaconda3/personal
source activate clustering

Rscript comparison_unweighted.R $PBS_ARRAY_INDEX

