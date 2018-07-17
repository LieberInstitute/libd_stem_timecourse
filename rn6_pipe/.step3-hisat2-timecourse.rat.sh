#!/bin/bash
#$ -cwd
#$ -l mem_free=5G,h_vmem=6G,h_fsize=200G
#$ -N step3-hisat2-timecourse.rat
#$ -pe local 8
#$ -o ./logs/hisat2-timecourse.$TASK_ID.txt
#$ -e ./logs/hisat2-timecourse.$TASK_ID.txt
#$ -t 1-36
#$ -tc 15
#$ -hold_jid pipeline_setup,step2-trim-timecourse.rat
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
echo "Sample id: $(cat /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")"
echo "****"

# Directories
mkdir -p /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/HISAT2_out/align_summaries
mkdir -p /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/HISAT2_out/infer_strandness

if [ FALSE == "TRUE" ]
then
    mkdir -p /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/HISAT2_out/unaligned
fi

## Locate file and ids
FILE1=$(awk 'BEGIN {FS="\t"} {print $1}' /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/samples.manifest | awk "NR==${SGE_TASK_ID}")
if [ TRUE == "TRUE" ] 
then
    FILE2=$(awk 'BEGIN {FS="\t"} {print $3}' /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/samples.manifest | awk "NR==${SGE_TASK_ID}")
fi
ID=$(cat /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")

if [ -f /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/trimmed_fq/${ID}_trimmed_forward_paired.fastq ] ; then
	## Trimmed, paired-end
	echo "HISAT2 alignment run on trimmed paired-end reads"
	FP=/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/trimmed_fq/${ID}_trimmed_forward_paired.fastq
	FU=/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/trimmed_fq/${ID}_trimmed_forward_unpaired.fastq
	RP=/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/trimmed_fq/${ID}_trimmed_reverse_paired.fastq
	RU=/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/trimmed_fq/${ID}_trimmed_reverse_unpaired.fastq
	
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/hisat2-2.0.4/hisat2 -p 8 	-x /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/ensembl/Rnor_6.0/hisat2_Rnor6.0toplevel -1 $FP -2 $RP -U ${FU},${RU} 	-S /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/HISAT2_out/${ID}_hisat_out.sam --rna-strandness RF --phred33      	2>/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/HISAT2_out/align_summaries/${ID}_summary.txt
	
elif  [ -f /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/trimmed_fq/${ID}_trimmed.fastq ] ; then
	## Trimmed, single-end
	echo "HISAT2 alignment run on trimmed single-end reads"
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/hisat2-2.0.4/hisat2 -p 8 	-x /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/ensembl/Rnor_6.0/hisat2_Rnor6.0toplevel -U /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/trimmed_fq/${ID}_trimmed.fastq 	-S /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/HISAT2_out/${ID}_hisat_out.sam --rna-strandness RF --phred33 	2>/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/HISAT2_out/align_summaries/${ID}_summary.txt

elif [ TRUE == "TRUE" ] ; then
	## Untrimmed, pair-end
	echo "HISAT2 alignment run on original untrimmed paired-end reads"
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/hisat2-2.0.4/hisat2 -p 8 	-x /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/ensembl/Rnor_6.0/hisat2_Rnor6.0toplevel -1 ${FILE1} -2 ${FILE2} 	-S /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/HISAT2_out/${ID}_hisat_out.sam --rna-strandness RF --phred33      	2>/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/HISAT2_out/align_summaries/${ID}_summary.txt

else
	## Untrimmed, single-end
	echo "HISAT2 alignment run on original untrimmed single-end reads"
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/hisat2-2.0.4/hisat2 -p 8 	-x /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/ensembl/Rnor_6.0/hisat2_Rnor6.0toplevel -U ${FILE1} 	-S /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/HISAT2_out/${ID}_hisat_out.sam --rna-strandness RF --phred33 	2>/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/HISAT2_out/align_summaries/${ID}_summary.txt
fi


###sam to bam
SAM=/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/HISAT2_out/${ID}_hisat_out.sam
ORIGINALBAM=/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/HISAT2_out/${ID}_accepted_hits.bam
SORTEDBAM=/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/HISAT2_out/${ID}_accepted_hits.sorted

#filter unmapped segments
echo "**** Filtering unmapped segments ****"
date
/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/samtools-1.2/samtools view -bh -F 4 ${SAM} > ${ORIGINALBAM}
/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/samtools-1.2/samtools sort -@ 8 ${ORIGINALBAM} ${SORTEDBAM}
/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/samtools-1.2/samtools index ${SORTEDBAM}.bam

## Clean up
rm ${SAM}
rm ${ORIGINALBAM}

## Run infer experiment
echo "**** Inferring strandedness with infer_experiment.py ****"
date
module load python/2.7.9
~/.local/bin/infer_experiment.py -i ${SORTEDBAM}.bam -r /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/RSeQC/rn6.bed 1> /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/HISAT2_out/infer_strandness/${ID}.txt 2>&1

echo "**** Job ends ****"
date
