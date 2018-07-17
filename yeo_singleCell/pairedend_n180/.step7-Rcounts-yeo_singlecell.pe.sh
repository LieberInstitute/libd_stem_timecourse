#!/bin/bash
#$ -cwd
#$ -pe local 5
#$ -l mem_free=18G,h_vmem=20G,h_fsize=200G
#$ -N step7-Rcounts-yeo_singlecell.pe
#$ -o ./logs/Rcounts-yeo_singlecell.txt
#$ -e ./logs/Rcounts-yeo_singlecell.txt
#$ -hold_jid pipeline_setup,step4-featCounts-yeo_singlecell.pe,step6-txQuant-yeo_singlecell.pe
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

Rscript /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/yeo_singleCell/pairedend_n180/.step7-create_count_objects-human.R -o hg38 -m /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/yeo_singleCell/pairedend_n180 -e yeo_singlecell -p pe -l TRUE -c FALSE -t 5 -s FALSE

echo "**** Job ends ****"
date
