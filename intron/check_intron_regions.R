# intron expression in human stem cells
library(derfinder)
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)

ensembl_v75 = import(con = "/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf", format = "gtf")
introns = import(con = '/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/Introns/introns_GRCh38.gtf', format = "gtf")

