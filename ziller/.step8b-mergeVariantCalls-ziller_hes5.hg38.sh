#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=8G,h_vmem=10G
#$ -N step8b-mergeVariantCalls-ziller_hes5.hg38
#$ -o ./logs/mergeVariantCalls-ziller_hes5.txt
#$ -e ./logs/mergeVariantCalls-ziller_hes5.txt
#$ -hold_jid pipeline_setup,step8-callVariants-ziller_hes5.hg38
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

VCFS=$(cat /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/ziller/samples.manifest | awk '{print "/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/ziller/Genotypes/"$NF".vcf.gz"}' | paste -sd " ")

module load vcftools
module load htslib
vcf-merge ${VCFS} | bgzip -c > /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/ziller/Genotypes/mergedVariants.vcf.gz

echo "**** Job ends ****"
date
