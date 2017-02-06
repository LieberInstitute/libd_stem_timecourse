
## load packages
library(jaffelab)
library(GenomicRanges)
library(rtracklayer)

#######################
# load genomic data ###
load("annotated_phenotype_data_stemCellTimecourse.rda")
load("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/rpkmCounts_AZpilot_jan6.hg38_n506.rda")
load("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/rawCounts_AZpilot_jan6.hg38_n506.rda")

## filter to same samples
geneCounts = geneCounts[,pd$SAMPLE_ID]
exonCounts = exonCounts[,pd$SAMPLE_ID]
jCounts = as.matrix(as.data.frame(jCounts[,gsub("-", ".", pd$SAMPLE_ID)]))

## filter junctions for expression levels
jIndex = which(rowSums(jCounts > 0) > 4)
jCounts = jCounts[jIndex,]
jMap = jMap[jIndex]

## add gene type
gtf = import("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf")
geneMap$gene_type = gtf$gene_type[match(geneMap$gencodeID, gtf$gene_id)]
exonMap$gene_type = gtf$gene_type[match(exonMap$gencodeID, gtf$gene_id)]
jMap$gene_type = gtf$gene_type[match(jMap$newGeneID, gtf$gene_id)]

## update junction annotation class
tt = rowSums(!is.na(as.data.frame(mcols(jMap)[,
	c("startExon", "endExon")])))
jMap$Class[tt == 1] = "AltStartEnd"
jMap$Class[tt == 2 & !jMap$inGencode] = "ExonSkip"

## add transcripts
txPath = "/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/Salmon_tx/"
txFiles = paste0(txPath, pd$SAMPLE_ID, "/quant.sf")
names(txFiles) = rownames(pd)
txList = lapply(txFiles, read.delim, row.names=1, as.is=TRUE)
tpmMat = sapply(txList, "[[", "TPM")
rownames(tpmMat) = ss(rownames(txList[[1]]),"|", fixed=TRUE)

############
## make expression sets
library(SummarizedExperiment)
pdDF = DataFrame(pd)
identical(pdDF$SAMPLE_ID, colnames(geneCounts))

## gene
geneMapGR = makeGRangesFromDataFrame(geneMap, keep=TRUE)
geneMapGR$gencodeTx = CharacterList(strsplit(geneMapGR$gencodeTx, ";"))
colnames(geneCounts) = rownames(pdDF)
rse_gene = SummarizedExperiment(
	assays = list('counts' = geneCounts),
    colData = pdDF, rowRanges = geneMapGR)

## gene
exonMapGR = makeGRangesFromDataFrame(exonMap, keep=TRUE)
exonMapGR$gencodeTx = CharacterList(strsplit(exonMapGR$gencodeTx, ";"))
colnames(exonCounts) = rownames(pdDF)
rse_exon = SummarizedExperiment(
	assays = list('counts' = exonCounts),
    colData = pdDF, rowRanges = exonMapGR)
	
## junction
jIndex = which(rowSums(jCounts > 0) > 2)
jCounts = jCounts[jIndex,]
jMap = jMap[jIndex]
colnames(jCounts) = rownames(pdDF)
rse_jxn = SummarizedExperiment(
	assays = list('counts' = jCounts),
    colData = pdDF, rowRanges = jMap)
	
## transcript
tx = gtf[which(gtf$type == "transcript")]
names(tx) = tx$transcript_id
txMap = tx[rownames(tpmMat)]
rse_tx = SummarizedExperiment(
	assays = list('tpm' = tpmMat),
    colData = pdDF, rowRanges = txMap)
	
#### function for RPKM
getRPKM = function(rse) {
	require(SummarizedExperiment)
	bg = matrix(rep(colData(rse)$totalMapped), 
		nc = ncol(rse), nr = nrow(rse),	byrow=TRUE)
	wid = matrix(rep(rowData(rse)$Length), 
		nr = nrow(rse), nc = ncol(rse),	byrow=FALSE)
	assays(rse)$counts/(wid/1000)/(bg/1e6)
}

getRPM = function(rse, target = 80e6) {
	require(SummarizedExperiment)
	bg = matrix(rep(colData(rse)$totalMapped/target), 
		nc = ncol(rse), nr = nrow(rse),	byrow=TRUE)
	assays(rse)$counts/bg
}

################
## save ########

save(rse_gene, getRPKM, compress=TRUE,
	file = "../data/libd_stemcell_timecourse_rseGene_n157.rda")
save(rse_exon, getRPKM, compress=TRUE,
	file = "../data/libd_stemcell_timecourse_rseExon_n157.rda")
save(rse_jxn, getRPM, compress=TRUE,
	file = "../data/libd_stemcell_timecourse_rseJxn_n157.rda")
save(rse_tx, compress=TRUE,
	file = "../data/libd_stemcell_timecourse_rseTx_n157.rda")
