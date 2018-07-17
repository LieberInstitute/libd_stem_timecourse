#!/bin/bash
#$ -cwd
#$ -N step3b-infer-strandness-scorecard.hg38
#$ -o ./logs/infer-strandness-scorecard.txt
#$ -e ./logs/infer-strandness-scorecard.txt
#$ -hold_jid pipeline_setup,step3-hisat2-scorecard.hg38
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

## Process the infer experiment info
Rscript /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/step3b_infer_strandness.R -o "HISAT2_out/infer_strandness" -p "inferred_strandness_pattern.txt"

echo "**** Job ends ****"
date
