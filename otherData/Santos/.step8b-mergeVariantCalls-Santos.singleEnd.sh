#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=8G,h_vmem=10G
#$ -N step8b-mergeVariantCalls-Santos.singleEnd
#$ -o ./logs/mergeVariantCalls-Santos.txt
#$ -e ./logs/mergeVariantCalls-Santos.txt
#$ -hold_jid pipeline_setup,step8-callVariants-Santos.singleEnd
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

VCFS=$(cat /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/otherData/Santos/samples.manifest | awk '{print "/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/otherData/Santos/Genotypes/"$NF".vcf.gz"}' | paste -sd " ")

module load vcftools
module load htslib
vcf-merge ${VCFS} | bgzip -c > /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/otherData/Santos/Genotypes/mergedVariants.vcf.gz

echo "**** Job ends ****"
date
