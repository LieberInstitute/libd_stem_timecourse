#!/bin/bash
#$ -cwd
#$ -pe local 5
#$ -l mem_free=18G,h_vmem=20G,h_fsize=200G
#$ -N step7-Rcounts-timecourse.rat
#$ -o ./logs/Rcounts-timecourse.txt
#$ -e ./logs/Rcounts-timecourse.txt
#$ -hold_jid pipeline_setup,step4-featCounts-timecourse.rat,step6-txQuant-timecourse.rat
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

Rscript /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/.step7-create_count_objects-rat.R -o rn6 -m /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe -e timecourse -p rat -l TRUE -c FALSE -t 5 -s reverse

echo "**** Job ends ****"
date
