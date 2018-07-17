#!/bin/bash
#$ -cwd
#$ -N step4-featCounts-Bardy.singleCell_clean
#$ -o ./logs/featCounts-Bardy_clean.txt
#$ -e ./logs/featCounts-Bardy_clean.txt
#$ -hold_jid pipeline_setup,step4-featCounts-Bardy.singleCell
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
rm -rf /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/otherData/Bardy/Counts/junction/tmpdir

echo "**** Job ends ****"
date
