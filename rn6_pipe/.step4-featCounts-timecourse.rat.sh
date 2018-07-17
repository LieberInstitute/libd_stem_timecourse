#!/bin/bash
#$ -cwd
#$ -l mem_free=8G,h_vmem=8G,h_fsize=100G
#$ -N step4-featCounts-timecourse.rat
#$ -pe local 8
#$ -o ./logs/featCounts-timecourse.$TASK_ID.txt
#$ -e ./logs/featCounts-timecourse.$TASK_ID.txt
#$ -t 1-36
#$ -tc 30
#$ -hold_jid pipeline_setup,step3-hisat2-timecourse.rat
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
mkdir -p /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/Counts/gene
mkdir -p /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/Counts/exon
mkdir -p /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/Counts/junction/tmpdir

FILE1=$(awk 'BEGIN {FS="\t"} {print $1}' /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/samples.manifest | awk "NR==${SGE_TASK_ID}")
if [ TRUE == "TRUE" ] 
then
    FILE2=$(awk 'BEGIN {FS="\t"} {print $3}' /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/samples.manifest | awk "NR==${SGE_TASK_ID}")
fi
ID=$(cat /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")
BAM=/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/HISAT2_out/${ID}_accepted_hits.sorted.bam

if [ TRUE == "TRUE" ] ; then 
	# genes	
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/subread-1.5.0-p3-source/bin/featureCounts 	-s 2 -p -T 8 -a /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/ensembl/Rnor_6.0/Rattus_norvegicus.Rnor_6.0.86.gtf 	-o /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/Counts/gene/${ID}_Ensembl.rnor6.0.rn6_Genes.counts $BAM
	# exons	
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/subread-1.5.0-p3-source/bin/featureCounts 	-s 2 -p -O -f -T 8 -a /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/ensembl/Rnor_6.0/Rattus_norvegicus.Rnor_6.0.86.gtf 	-o /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/Counts/exon/${ID}_Ensembl.rnor6.0.rn6_Exons.counts $BAM
else
	# genes	
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/subread-1.5.0-p3-source/bin/featureCounts 	-s 2 -T 8 -a /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/ensembl/Rnor_6.0/Rattus_norvegicus.Rnor_6.0.86.gtf 	-o /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/Counts/gene/${ID}_Ensembl.rnor6.0.rn6_Genes.counts $BAM
	# exons	
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/subread-1.5.0-p3-source/bin/featureCounts 	-s 2 -O -f -T 8 -a /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/ensembl/Rnor_6.0/Rattus_norvegicus.Rnor_6.0.86.gtf 	-o /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/Counts/exon/${ID}_Ensembl.rnor6.0.rn6_Exons.counts $BAM
fi
	
# junctions	
OUTJXN=/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/Counts/junction/${ID}_junctions_primaryOnly_regtools.bed
OUTCOUNT=/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/Counts/junction/${ID}_junctions_primaryOnly_regtools.count
TMPDIR=/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/Counts/junction/tmpdir
TMPBAM=${TMPDIR}/${ID}.bam
#filter only primary alignments
/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/samtools-1.2/samtools view -@ 8 -bh -F 0x100 $BAM > ${TMPBAM}
/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/samtools-1.2/samtools index ${TMPBAM}

## Load python 2.7.9 since the default one cannot run:
# python
# import site
module load python/2.7.9
/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/regtools/build/regtools junctions extract -i 9 -o ${OUTJXN} ${TMPBAM}
/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/bed_to_juncs_withCount < ${OUTJXN} > ${OUTCOUNT}


echo "**** Job ends ****"
date
