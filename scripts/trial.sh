#!/bin/bash

#SBATCH --job-name='hlatyping'
#SBATCH --output=hlatyping-%j-stdout.log
#SBATCH --error=hlatyping-%j-stderr.log
#SBATCH --time=7-00:00:00
#SBATCH --mem=80GB
##SBATCH --partition=HighMem
##SBATCH --nodelist=highmem-002
##SBATCH --nodes=20
##SBATCH --cpus-per-task=20

echo "Submitting SLURM job"
nextflow run main.nf -profile slurm -resume 
