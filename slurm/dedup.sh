#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --time=1-00:00:00

file=/projects/bgmp/czakari2/bioinformatics/Bi624/Deduper-czakarian/C1_SE_uniqAlign_sorted.sam
umi=/projects/bgmp/czakari2/bioinformatics/Bi624/Deduper-czakarian/STL96.txt
script=/projects/bgmp/czakari2/bioinformatics/Bi624/Deduper-czakarian/zakarian_deduper.py


/usr/bin/time -v $script -f $file -u $umi
