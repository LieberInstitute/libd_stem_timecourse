#Use Gencode v25 GTF (/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf) to get strand-specific counts of all introns [maybe unique of all `GenomicFeatures::intronsByTranscript()`] plus library-size normalized bigwigs (/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/Coverage/) for N=139 - want to show loss of intronic expression over differentiation
#Badoi and Steve Stem-cell 

#`/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/junction_txdb/make_txdb_from_gtf.R`
#We can make a modified form of that to create an intron GTF file to input into FeatureCounts
#######

library(GenomicFeatures)
library(GenomicRanges)
library(org.Hs.eg.db)
library(rtracklayer)
source("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/eqtl_functions.R")

################################################################
####  hg19  ####################################################
################################################################

## chromosome info
chrInfo = getChromInfoFromUCSC("hg19")
chrInfo$chrom = as.character(chrInfo$chrom)
chrInfo = chrInfo[chrInfo$chrom %in% paste0("chr", c(1:22,"X","Y", "M")),]
chrInfo$isCircular = rep(c(FALSE, TRUE), c(24,1))
si = with(chrInfo, Seqinfo(as.character(chrom), length, isCircular, genome="hg19"))

##########################
##### ENSEMBL ############

## read in GTF as GRanges
ensembl_v75 = import(con = "/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf", format = "gtf")
#seqlevels(ensembl_v75,force=TRUE) = c(1:22,"X","Y","MT")
#seqlevels(ensembl_v75) = paste0("chr", c(1:22,"X","Y","M"))
#seqinfo(ensembl_v75) = si

# get map
ensMap = mcols(ensembl_v75)
ensMap = ensMap[!duplicated(ensMap$transcript_id),]

##  convert to txdb
ensembl_v75_txdb = makeTxDbFromGRanges(ensembl_v75)
#saveDb(ensembl_v75_txdb, file="ensembl_v75.sqlite")

# get introns
introns = intronsByTranscript(ensembl_v75_txdb,use.names=TRUE)
introns = unlist(introns)
introns$TranscriptID = names(introns)
#introns$ensemblID = ensMap$gene_id[match(introns$TranscriptID, ensMap$transcript_id)]
#introns$symbol = ensMap$gene_name[match(introns$TranscriptID, ensMap$transcript_id)]
mcols(introns) <-ensMap[match(names(introns), ensMap$transcript_id),]

unique_introns <- unique(introns)
export(unique_introns, '/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/Introns/introns_GRCh38.gtf', format = 'gtf')
