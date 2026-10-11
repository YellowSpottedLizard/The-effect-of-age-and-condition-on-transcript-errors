#!/bin/bash
#SBATCH -A naiss2025-5-81
#SBATCH -p main
#SBATCH -t 5:00:00
#SBATCH -J getregion

set -e

module load PDC
module load freebayes
module load samtools

# bam file
sample=$1 # ./bam/P20702_102.Aligned.sortedByCoord.out.bam
sampleID=$(awk -F "." '{ print $1 }' <<< $(basename $sample)) # P20702_102

refindex=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/refgenome/genome.fa.fai
samtools depth $sample | coverage_to_regions.py ${refindex} 300 > $sample.region
