
## load packages
library(jaffelab)
library(GenomicRanges)
library(rtracklayer)

ERCCDIR = '/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/Ercc'

## get pd for subset of samples, n=157
load("../prep_samples/annotated_phenotype_data_stemCellTimecourse.rda")
sampIDs = as.vector(pd$SAMPLE_ID)

## observed kallisto tpm
erccMat = sapply(sampIDs, function(x) {
  read.table(file.path(ERCCDIR, x, "abundance.tsv"),header = TRUE)$tpm
})
rownames(erccMat) = read.table(file.path(ERCCDIR,sampIDs[1],"abundance.tsv"),header = TRUE)$target_id
colnames(erccMat) = pd$SampleID

erccLength = read.table(file.path(ERCCDIR,sampIDs[1],"abundance.tsv"),header = TRUE)$length
erccLength = data.frame(ercc=as.character(rownames(erccMat)),length=erccLength)
erccLength$ercc = as.character(erccLength$ercc)

############
## make expression sets
library(SummarizedExperiment)
pdDF = DataFrame(pd)

## make rse
rse_ercc = SummarizedExperiment(
	assays = list('tpm' = erccMat),
    colData = pdDF, rowData=erccLength)

						
################
## save ########

save(rse_ercc, compress=TRUE,
	file = "libd_stemcell_timecourse_rseERCC_n157.rda")
