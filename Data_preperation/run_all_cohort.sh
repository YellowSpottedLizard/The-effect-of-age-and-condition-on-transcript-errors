#!/bin/bash
#SBATCH -A naiss2025-5-81
#SBATCH -p shared
#SBATCH -t 48:00:00
#SBATCH -J run_all_cohort
set -e

codedir=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/LTS-U_Friberg_2205-buffering/code/TE
bam_dir=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/bam
output1=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/outputs
output2=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/outputs_rackham

first=1
for file in ${bam_dir}/P20702_*.Aligned.sortedByCoord.out.bam.region; do 
	sampleID=$(awk -F "." '{ print $1 }' <<< $(basename $file))
	echo ${sampleID}
	#if [ ! -f outputs_rackham/${sampleID}_Y* ]; then
	if [ ! -f ${output1}/${sampleID}_Y* ] && [ ! -f ${output2}/${sampleID}_Y* ]; then
		if [ ${first} == 1 ]; then
				call=$(sbatch ${codedir}/run_all.sh $file)
				jobid=$(awk '{ print $4 }' <<< $call)
				echo "submitted ${sampleID} with jobid ${jobid}"
			first=0
		else 
			echo "waiting for ${jobid}"
			slurmfile=slurm-${jobid}.out
			until [ -f slurm-${jobid}.out ]
			do
			    sleep 30
			done
			echo "job ${jobid} started"	
			region_jobids=$(<slurm-${jobid}.out)
			echo "Freebayes region jobids: ${region_jobids}"	
			call=$(sbatch --dependency=afterok:${region_jobids} ${codedir}/run_all.sh $file)
                                jobid=$(awk '{ print $4 }' <<< $call)
				echo "submitted ${jobid}"
		fi
	else
	    echo "skipping ${sampleID}"
	fi

done
