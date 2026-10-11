#!/bin/bash
#SBATCH -A naiss2025-5-81
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH -J run_all


set -e

# bam file
regions=$1 # ./bam/P20702_102.Aligned.sortedByCoord.out.bam.region
sampleID=$(awk -F "." '{ print $1 }' <<< $(basename $regions)) # P20702_102
codedir=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/LTS-U_Friberg_2205-buffering/code/TE

jobids=''

while read region; do
	# 144 GB for 12 cores
	#sbatch ./run_freebayes_region.sh $sampleID $region
	call=$(sbatch ${codedir}/run_freebayes_region.sh $sampleID $region)
	jobid=$(awk '{ print $4 }' <<< $call)
	jobids="${jobids}:${jobid}"
done<$regions

echo ${jobids:1}
