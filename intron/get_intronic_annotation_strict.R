#Variation on code from: https://github.com/LieberInstitute/RNAseq-pipeline/commit/323e6b1405322b5674a21800f2af969c99d3aeef
library('GenomicRanges')
library('GenomicFeatures')
library('bumphunter')
#Getting updated functions for annotation written by Leo
source("./annotateTranscripts_fix.R")
source("./matchGenes_fix.R")


load('/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/gs/gs_gencode_v25_hg38.Rdata')
strict_introns <- unique( gs_gencode_v25_hg38$fullGenome[gs_gencode_v25_hg38$fullGenome$theRegion == 'intron'] )

### Annotating strict_introns to nearest gene

## Get the chromosome info for hg38
 chrInfo <- getChromInfoFromUCSC('hg38')
 chrInfo$chrom <- as.character(chrInfo$chrom)
 chrInfo <- chrInfo[chrInfo$chrom %in% paste0('chr', c(1:22, 'X', 'Y', 'M')), ]
 chrInfo$isCircular <- rep(c(FALSE, TRUE), c(24, 1))
 si <- with(chrInfo, Seqinfo(as.character(chrom), length, isCircular,
     genome = 'hg38'))
 
 
###############################
### Import Gencode Info #######
gencode_v25 <- rtracklayer::import('/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf')

## Add the chromosome lengths
 seqlevels(gencode_v25,force=TRUE) <- paste0('chr', c(1:22, 'X', 'Y', 'M'))
 seqinfo(gencode_v25) <- si
gencode_v25_txdb <- makeTxDbFromGRanges(gencode_v25)

genes <- annotateTranscripts_fix(gencode_v25_txdb,annotationPackage = 'org.Hs.eg.db' )

ann <- matchGenes_fix(as.data.frame(strict_introns), subject = genes)
mcols(strict_introns) = cbind(mcols(strict_introns), ann)

#supportedUCSCtables(genome="hg38")

save(strict_introns,file='/users/ssemick/StemCell/rdas/strict_introns_annotation.rda')
table(strict_introns$subregion)
strict_introns$strand =NULL
#rtracklayer::export(strict_introns, '/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/Introns/introns_GRCh38_annotated.gtf', format = 'gtf')
#rtracklayer::export(strict_introns, '/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/Introns/introns_GRCh38_annotated.bed', format = 'bed')

#
test = sample(strict_introns,100)
start(strict_introns) = start(test)+1
end(test) = end(test)-1
ann <- matchGenes_fix(as.data.frame(test), subject = genes)
table(ann$subregion)
test = test[which(ann$subregion=='covers exon(s)')]
mcols(test) = mcols(test)$theRegion

rtracklayer::export(test, '/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/Introns/test_introns.gtf', format = 'gtf')

