#!/bin/bash
#$ -cwd
#$ -pe local 5
#$ -l bluejay,mem_free=18G,h_vmem=20G,h_fsize=200G
#$ -N step7-Rcounts-ziller_hes5.hg38
#$ -o ./logs/Rcounts-ziller_hes5.txt
#$ -e ./logs/Rcounts-ziller_hes5.txt
#$ -hold_jid pipeline_setup,step4-featCounts-ziller_hes5.hg38,step6-txQuant-ziller_hes5.hg38
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

Rscript /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/ziller/.step7-create_count_objects-human.R -o hg38 -m /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/ziller -e ziller_hes5 -p hg38 -l TRUE -c FALSE -t 5 -s FALSE

echo "**** Job ends ****"
date
