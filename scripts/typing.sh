#!/bin/bash

## SBATCH --job-name='hlatyping'
## SBATCH --cpus-per-task=20
## SBATCH --output=hlatyping-%j-stdout.log
## SBATCH --error=hlatyping-%j-stderr.log
##SBATCH --time=14-00:00:00
##SBATCH --partition=HighMem
##SBATCH --nodelist=highmem-002
##SBATCH --nodes=20

echo "Submitting SLURM job"
for id in `cat h3a.txt`
do
#echo Now running ${id}
HLA-LA.pl --BAM /cbio/projects/013/custom-bam.ruth/selected/${id}/${id}*.bam \
--graph PRG_MHC_GRCh38_withIMGT --sampleID ${id} \
--workingDir ./output --maxThreads 20 &
done

