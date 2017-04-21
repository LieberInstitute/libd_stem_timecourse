INTRON=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/Introns/test_introns.gtf
EXON=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf

#module load bedtools

#determine if there is 1 base overlap of intron gtf w/ exon gtf
bedtools intersect -wao -a $INTRON -b $EXON  > testIntronOverlap.bed