#!/bin/bash
#$ -cwd
#$ -pe local 5
#$ -l bluejay,mem_free=18G,h_vmem=20G,h_fsize=200G
#$ -N step7-Rcounts-Bardy.singleCell
#$ -o ./logs/Rcounts-Bardy.txt
#$ -e ./logs/Rcounts-Bardy.txt
#$ -hold_jid pipeline_setup,step4-featCounts-Bardy.singleCell,step6-txQuant-Bardy.singleCell
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

## Force R 3.3.x in JHPCE (to avoid some issues with conda_R)
module unload conda_R
module load R/3.3.x

Rscript /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/otherData/Bardy/.step7-create_count_objects-human.R -o hg38 -m /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/otherData/Bardy -e Bardy -p singleCell -l TRUE -c TRUE -t 5 -s FALSE

echo "**** Job ends ****"
date
