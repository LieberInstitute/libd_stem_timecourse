#!/bin/bash
#$ -cwd
#$ -l mem_free=8G,h_vmem=10G
#$ -N step8b-mergeVariantCalls-yeo_singlecell.pe
#$ -o ./logs/mergeVariantCalls-yeo_singlecell.txt
#$ -e ./logs/mergeVariantCalls-yeo_singlecell.txt
#$ -hold_jid pipeline_setup,step8-callVariants-yeo_singlecell.pe
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

VCFS=$(cat /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/yeo_singleCell/pairedend_n180/samples.manifest | awk '{print "/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/yeo_singleCell/pairedend_n180/Genotypes/"$NF".vcf.gz"}' | paste -sd " ")

module load vcftools
module load htslib
vcf-merge ${VCFS} | bgzip -c > /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/yeo_singleCell/pairedend_n180/Genotypes/mergedVariants.vcf.gz

echo "**** Job ends ****"
date
