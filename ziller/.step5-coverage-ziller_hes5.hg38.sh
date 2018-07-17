#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=12G,h_vmem=16G,h_fsize=100G
#$ -N step5-coverage-ziller_hes5.hg38
#$ -o ./logs/coverage-ziller_hes5.$TASK_ID.txt
#$ -e ./logs/coverage-ziller_hes5.$TASK_ID.txt
#$ -t 1-19
#$ -tc 40
#$ -hold_jid pipeline_setup,step3-hisat2-ziller_hes5.hg38,step3b-infer-strandness-ziller_hes5.hg38
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"
echo "Sample id: $(cat /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/ziller/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")"
echo "****"

if [ ! -f "inferred_strandness_pattern.txt" ]
then
    echo "Missing the file inferred_strandness_pattern.txt"
    exit 1
fi

FILE1=$(awk 'BEGIN {FS="\t"} {print $1}' /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/ziller/samples.manifest | awk "NR==${SGE_TASK_ID}")
if [[ TRUE == "TRUE" ]]
then
    FILE2=$(awk 'BEGIN {FS="\t"} {print $3}' /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/ziller/samples.manifest | awk "NR==${SGE_TASK_ID}")
fi
ID=$(cat /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/ziller/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")
STRANDRULE=$(cat inferred_strandness_pattern.txt)

## Can only use -d when the data is stranded
if [ ${STRANDRULE} == "none" ]
then
    STRANDPARAM=""
    STRANDOPTION=""
else
    STRANDPARAM="-d "
    STRANDOPTION="\"${STRANDRULE}\""
fi

## Normalizing bigwigs to 40 million 100 bp reads
module load python/2.7.9
module load ucsctools
python ~/.local/bin/bam2wig.py -s /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg38.chrom.sizes.gencode -i /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/ziller/HISAT2_out/${ID}_accepted_hits.sorted.bam -t 4000000000 -o /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/ziller/Coverage/${ID}  

## Remove temp files
rm /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/ziller/Coverage/${ID}*.wig

echo "**** Job ends ****"
date
