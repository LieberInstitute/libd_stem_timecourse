#!/bin/bash
#$ -cwd
#$ -pe local 5
#$ -l bluejay,mem_free=18G,h_vmem=20G,h_fsize=200G
#$ -N step7-Rcounts-scorecard.hg38
#$ -o ./logs/Rcounts-scorecard.txt
#$ -e ./logs/Rcounts-scorecard.txt
#$ -hold_jid pipeline_setup,step4-featCounts-scorecard.hg38,step6-txQuant-scorecard.hg38
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

Rscript /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/scorecard/.step7-create_count_objects-human.R -o hg38 -m /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/scorecard -e scorecard -p hg38 -l TRUE -c FALSE -t 5 -s FALSE

echo "**** Job ends ****"
date
