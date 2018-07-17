
## load packages
library(jaffelab)
library(GenomicRanges)
library(SummarizedExperiment)
library(rtracklayer)

#######################
# load genomic data ###
load("annotated_phenotype_data_stemCellTimecourse.rda")

load("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/rse_gene_AZpilot_jan6.hg38_n506.Rdata")
geneCounts = assays(rse_gene)$counts
geneMap = as.data.frame(rowRanges(rse_gene))
names(geneMap)[1:5] = c("Chr","Start","End","Width","Strand")
rm(rse_gene)

load("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/rse_exon_AZpilot_jan6.hg38_n506.Rdata")
exonCounts = assays(rse_exon)$counts
exonMap = as.data.frame(rowRanges(rse_exon))
names(exonMap)[1:5] = c("Chr","Start","End","Width","Strand")
rm(rse_exon)

load("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/rse_jx_AZpilot_jan6.hg38_n506.Rdata")
jCounts = as.data.frame(assays(rse_jx)$counts)
jMap = rowRanges(rse_jx)
rm(rse_jx)


## filter to 146 samples
dropIndex = which(pd$LINE == "165-B-8X" | pd$RNA_NO %in% c("R16-081", "R16-094", "R16-095", "R16-197", "R16-288"))
pd = pd[-dropIndex,]
geneCounts = geneCounts[,pd$SAMPLE_ID]
exonCounts = exonCounts[,pd$SAMPLE_ID]
jCounts = as.matrix(jCounts[,gsub("-", ".", pd$SAMPLE_ID)])

## filter junctions for expression levels
jIndex = which(rowSums(jCounts > 0) > 4)
jCounts = jCounts[jIndex,]
jMap = jMap[jIndex]


# ## add gene type
# gtf = import("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf")
# geneMap$gene_type = gtf$gene_type[match(geneMap$gencodeID, gtf$gene_id)]
# exonMap$gene_type = gtf$gene_type[match(exonMap$gencodeID, gtf$gene_id)]
# jMap$gene_type = gtf$gene_type[match(jMap$newGeneID, gtf$gene_id)]

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

## exon
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


## add mean expression
jRpkm = getRPM(rse_jxn)
jMap$meanExprs = rowMeans(jRpkm)
rse_jxn = SummarizedExperiment(
	assays = list('counts' = jCounts),
    colData = pdDF, rowRanges = jMap)

	
################
## save ########

save(rse_gene, getRPKM, compress=TRUE,
	file = "../data/libd_stemcell_timecourse_rseGene_n146.rda")
save(rse_exon, getRPKM, compress=TRUE,
	file = "../data/libd_stemcell_timecourse_rseExon_n146.rda")
save(rse_jxn, getRPM, compress=TRUE,
	file = "../data/libd_stemcell_timecourse_rseJxn_n146.rda")
save(rse_tx, compress=TRUE,
	file = "../data/libd_stemcell_timecourse_rseTx_n146.rda")
