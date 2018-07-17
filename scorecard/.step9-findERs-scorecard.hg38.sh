#!/bin/bash
#$ -cwd
#$ -pe local 10
#$ -l bluejay,mem_free=15G,h_vmem=15G,h_fsize=200G
#$ -N step9-findERs-scorecard.hg38
#$ -o ./logs/findERs-scorecard.txt
#$ -e ./logs/findERs-scorecard.txt
#$ -hold_jid pipeline_setup,step5b-meanCoverage-scorecard.hg38
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

for meanFile in Coverage/mean*.bw
do
    echo "************************************"
    date
    echo "Initializing script for ${meanFile}"
    echo "************************************"
    Rscript /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/scorecard/.step9-find_expressed_regions.R -m ${meanFile} -o /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/scorecard/ERs -i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg38.chrom.sizes.gencode -c 10
done

echo "**** Job ends ****"
date
