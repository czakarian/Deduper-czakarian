#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --time=1-00:00:00

#file=/projects/bgmp/shared/deduper/C1_SE_uniqAlign.sam
file=/projects/bgmp/shared/deduper/C1_PE_uniqAligned.sam

#out=C1_SE_uniqAlign_sorted.sam
out=C1_PE_uniqAligned_sorted.sam

conda activate bgmp_py39

cd ..

/usr/bin/time -v samtools sort $file -o $out