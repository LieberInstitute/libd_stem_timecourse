#!/bin/bash
#$ -cwd
#$ -pe local 5
#$ -l bluejay,mem_free=18G,h_vmem=20G,h_fsize=200G
#$ -N step7-Rcounts-yeo_singlecell.se
#$ -o ./logs/Rcounts-yeo_singlecell.txt
#$ -e ./logs/Rcounts-yeo_singlecell.txt
#$ -hold_jid pipeline_setup,step4-featCounts-yeo_singlecell.se,step6-txQuant-yeo_singlecell.se
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

Rscript /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/yeo_singleCell/singleend_n34/.step7-create_count_objects-human.R -o hg38 -m /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/yeo_singleCell/singleend_n34 -e yeo_singlecell -p se -l FALSE -c FALSE -t 5 -s FALSE

echo "**** Job ends ****"
date
