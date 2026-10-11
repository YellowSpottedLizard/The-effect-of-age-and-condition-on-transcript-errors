#!/bin/bash
#SBATCH -A naiss2025-5-81
#SBATCH -p shared
#SBATCH -c 50
#SBATCH -t 5:00:00
#SBATCH -J freebayes

set -e

module load PDC
#module load bioinfo-tools
module load freebayes

sample=$1
region=$2
bam_dir=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/bam
ref=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/refgenome/genome.fa
output=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/outputs

echo ${sample} ${region}

freebayes -f ${ref} \
-r ${region} \
--haplotype-length 0 \
--min-alternate-count 1 \
--min-alternate-fraction 0 \
--pooled-continuous \
--min-coverage 100 \
--use-best-n-alleles 4 \
-i \
${bam_dir}/${sample}.Aligned.sortedByCoord.out.bam > ${output}/${sample}_${region}.freebayes.vcf

echo "###" >> ${output}/${sample}_${region}.freebayes.vcf
