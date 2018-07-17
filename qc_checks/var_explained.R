### load packages
library(jaffelab)
library(GenomicRanges)
library(SummarizedExperiment)
library(RColorBrewer)
library(rafalib)

MAINDIR = '/dcl01/lieber/ajaffe/lab/libd_stem_timecourse'
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n146.rda"))

geneCounts = assays(rse_gene)$counts
geneRpkm = getRPKM(rse_gene)
geneMap = rowData(rse_gene)

pd = colData(rse_gene)
pd$REP = as.numeric(pd$BIO_REP)
pd$REP[is.na(pd$REP)] = 1
pd$DAY = ordered(pd$DAY, levels=sort(as.numeric(levels(as.factor(pd$DAY)))) )
pd$DIV = as.numeric(as.character(pd$DAY))
pd$EXPERIMENT[pd$EXPERIMENT=="TIMECOURSE1"] = "TIMECOURSE"
pd$LINE = gsub('-','_',pd$LINE)

## don't use technical replicates
bioInd = which(pd$Class == "Naked genomes")
pd = pd[bioInd,]
geneRpkm = geneRpkm[,bioInd]
geneCounts = geneCounts[,bioInd]   # n=128

## don't use RENEW controls or NEURONS_ALONE in analyses
dropInd = which(pd$CONDITION %in% c("RENEW","NEURONS_ALONE"))
pd = pd[-dropInd,]
geneRpkm = geneRpkm[,-dropInd]
geneCounts = geneCounts[,-dropInd]   # n=106


pd$COND = pd$CONDITION
pd$COND[grep("NEURON",pd$COND)] = "NEURON"  ## rename NEURONS_PLUS_ASTROS
pd$COND[grep("ACC_DORSAL",pd$COND)] = "ACC_DORSAL"  ## rename ACC_DORSAL(2)
pd$COND = as.factor(pd$COND)
pd$COND = factor(pd$COND,levels(pd$COND)[c(1,3,4,2)])  ## put levels in order


###########################
#### filter expression ####
geneMap$meanExprs = rowMeans(geneRpkm)
gIndex = which(geneMap$meanExprs > 0.1)
geneRpkm = geneRpkm[gIndex,]
geneCounts = geneCounts[gIndex,]
geneMap = geneMap[gIndex,]
map = geneMap

## merge data
yExprs = as.matrix(log2(geneRpkm+1))

##### variance explained
library(doMC)
registerDoMC(cores=12)

varCompAnalysis = foreach(i = 1:nrow(yExprs)) %dopar% {
	if(i %% 1000 == 0) cat(".")
	fit = lm(yExprs[i,]	~ factor(DIV)  +DX +  Donor + LINE + 
		totalAssignedGene + ERCCsumLogErr, data=pd)
	full = anova(fit)
	fullSS =full$"Sum Sq"
	signif(cbind(full,PctExp=fullSS/
		sum(fullSS)*100),3)
}
names(varCompAnalysis) = rownames(yExprs)
varCompMat = t(sapply(varCompAnalysis, "[[", "PctExp"))
colnames(varCompMat) = rownames(varCompAnalysis[[1]])
boxplot(varCompMat)