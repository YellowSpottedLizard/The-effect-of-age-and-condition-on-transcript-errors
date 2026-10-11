#!/bin/bash
#SBATCH -A naiss2025-5-81
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -t 10:00:00
#SBATCH -J collect_freebayes_results
set -e

module load PDC
module load bedtools

#This script should be run in /cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors
#It checks if freebayes is complete for all regions for a sample, 
#And if so concatentes the freebayes region files and put them into a specific folder

codedir=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/LTS-U_Friberg_2205-buffering/code/TE
runlog=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/freebayes_runlog_20250804.txt
concat_dir=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/concatenated_freebayes
outputs=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/outputs
outputs_rackham=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/outputs_rackham
bam=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/bam
genes=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/refgenome/genes_dm6.bed
intergenic=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/refgenome/intergenic_dm6.bed

echo 'Freebayes runlog '$(date) > ${runlog}

for sample in bam/P20702_*.Aligned.sortedByCoord.out.bam.region; do 
	sampleID=$(awk -F "." '{ print $1 }' <<< $(basename $sample))
	echo ${sampleID}

	#First check if sample already concatenated
	if [ ! -f ${concat_dir}/${sampleID}.concatenated.vcf ]; then

		#check if results are in the outputs folder:
		if [ -f ${outputs}/${sampleID}_Y* ]; then
			echo "processed on Dardel"
			#Check how many regions we should expect
			n_regions=$(wc -l < ${bam}'/'${sampleID}'.Aligned.sortedByCoord.out.bam.region')
			echo "${n_regions}"
			#Check how many complete vcf files we have:
			n_vcfs=$(tail -n 1 ${outputs}'/'${sampleID}*'.freebayes.vcf' | grep '###' | wc -l)
			echo ${n_vcfs}

			if [ ${n_regions} == ${n_vcfs} ]; then
				echo ${sampleID}' '${n_regions}' '${n_vcfs}' Done' >> ${runlog} 
				#concatenate
				cat ${outputs}/${sampleID}*vcf >> ${concat_dir}'/'${sampleID}.concatenated.vcf
			else
				echo ${sampleID}' '${n_regions}' '${n_vcfs}' output Restart' >> ${runlog} 
			fi

		#Second checks if sample was processed on Rackham		
		elif [ -f ${outputs_rackham}/${sampleID}_Y* ]; then

			echo "processed on Rackham"
			#Check how many regions we should expect
			n_regions=$(wc -l < ${bam}'/'${sampleID}'.Aligned.sortedByCoord.out.bam.region')
			echo ${n_regions}
			#Check how many complete vcf files we have:
			n_vcfs=$(tail -n 1 outputs_rackham/${sampleID}*.freebayes.vcf | grep '###' | wc -l)
			echo ${n_vcfs}
			if [ ${n_regions} == ${n_vcfs} ]; then
				echo ${sampleID}' '${n_regions}' '${n_vcfs}' Done' >> ${runlog} 
				#concatenate
				cat ${outputs_rackham}/${sampleID}*vcf >> ${concat_dir}'/'${sampleID}.concatenated.vcf
			else
				echo ${sampleID}' '${n_regions}' '${n_vcfs}' output_rackham Restart' >> ${runlog} 
			fi

		else 
			echo ${sampleID}' '${n_regions}' '${n_vcfs}' Not_run' >> ${runlog} 
		fi
		#extracts variants in genes and in intergenic regions
		#Note -u is used to print variants only once even if a avriant overlaps several genes. This happens since genes can overlap each other. 
		bedtools intersect -a ${concat_dir}'/'${sampleID}.concatenated.vcf -b ${genes} -header -u > ${concat_dir}'/'${sampleID}.genes.vcf
		bedtools intersect -a ${concat_dir}'/'${sampleID}.concatenated.vcf -b ${intergenic} -header -u > ${concat_dir}'/'${sampleID}.intergenic.vcf
	
	fi
	
done
