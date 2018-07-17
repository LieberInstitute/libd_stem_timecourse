#!/bin/bash
#$ -cwd
#$ -N step4-featCounts-yeo_singlecell.pe_clean
#$ -o ./logs/featCounts-yeo_singlecell_clean.txt
#$ -e ./logs/featCounts-yeo_singlecell_clean.txt
#$ -hold_jid pipeline_setup,step4-featCounts-yeo_singlecell.pe
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

## Delete temporary files after they have been used
rm -rf /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/yeo_singleCell/pairedend_n180/Counts/junction/tmpdir

echo "**** Job ends ****"
date
