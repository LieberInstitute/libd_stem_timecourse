# !/bin/sh

# qsub -cwd -l mf=10G,h_vmem=18G -t 2-74 -tc 20 download_reads.sh
# qsub -cwd -pe local 3 -l mf=5G,h_vmem=8G -t 2-74 download_reads.sh

# awk -F"\t" '{print $6}' scorecard_SraRunTable.txt > srrs.txt
FILELIST=/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/scorecard/srrs.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)
fastq-dump -O /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/scorecard/FASTQ/ --gzip --split-files $ID 

# R=/dcl01/lieber/ajaffe/PublicData/BloodDegrade/FASTQ/${ID}.fastq.gz
# O=/dcl01/lieber/ajaffe/PublicData/BloodDegrade/TopHat/$ID
# tophat2 -p 3 -G /dcl01/lieber/ajaffe/Annotation/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf -o $O /dcl01/lieber/ajaffe/Annotation/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome $R
# BAM=/dcl01/lieber/ajaffe/PublicData/BloodDegrade/TopHat/$ID/accepted_hits.bam
# samtools index $BAM

# ## feature count
# featureCounts -A /dcl01/lieber/ajaffe/Annotation/chrAliases_GRCh37_to_hg19.csv \
	# -a /dcl01/lieber/ajaffe/Annotation/Homo_sapiens.GRCh37.75.gtf \
	# -o /dcl01/lieber/ajaffe/PublicData/BloodDegrade/Counts/Gene/${ID}_Ensembl_v75_Genes.counts $BAM
# featureCounts -O -f -A /dcl01/lieber/ajaffe/Annotation/chrAliases_GRCh37_to_hg19.csv \
	# -a /dcl01/lieber/ajaffe/Annotation/Homo_sapiens.GRCh37.75.gtf \
	# -o /dcl01/lieber/ajaffe/PublicData/BloodDegrade/Counts/Exon/${ID}_Ensembl_v75_Exons.counts $BAM

# # junctions	
# OUTJXN=/dcl01/lieber/ajaffe/PublicData/BloodDegrade/Junctions/${ID}_junctions_primaryOnly_regtools.bed
# TMPBAM=$TMPDIR/${ID}.bam
# samtools view -bh -F 0x100 $BAM > $TMPBAM
# samtools index $TMPBAM
# regtools junctions extract -i 9 -o $OUTJXN $TMPBAM
# OUTCOUNT=/dcl01/lieber/ajaffe/PublicData/BloodDegrade/Junctions/${ID}_junctions_primaryOnly_regtools.count
# /users/ajaffe/Lieber/Projects/RNAseq/bed_to_juncs_withCount < $OUTJXN > $OUTCOUNT

# bw
BG=/dcl01/lieber/ajaffe/PublicData/BloodDegrade/Coverage/${ID}.bedGraph
CHRSIZE=/users/ajaffe/Lieber/Projects/RNAseq/hg19.chrom.sizes
BW=/dcl01/lieber/ajaffe/PublicData/BloodDegrade/Coverage/${ID}.bw

bedtools genomecov -ibam $BAM -bga -split > $BG
awk '$1 ~ /^chr/' $BG > ${BG}.tmp
bedGraphToBigWig ${BG}.tmp $CHRSIZE $BW
rm $BG ${BG}.tmp
