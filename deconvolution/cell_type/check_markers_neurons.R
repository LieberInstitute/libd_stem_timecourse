######

library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(genefilter)
library(RColorBrewer)
library(limma)
library(edgeR)

pal = brewer.pal(10, "Set3")
pal[2] = "#8B795E" # darken NPCs
 
MAINDIR = "/dcl01/lieber/ajaffe/lab/libd_stem_timecourse"

# load gene-level
load("/dcl01/ajaffe/data/lab/singleCell/Darmanis/rna-seq-pipeline_run2/rse_gene_Darmanis_scRNASeq_Darmanis_n466.Rdata")
rse_geneQuake = rse_gene
rm(rse_gene)

# phenotype
load("/dcl01/lieber/ajaffe/PublicData/Brain/Darmanis_Quake/rdas/phenotype_quake_processed.rda")
colData(rse_geneQuake) = cbind(colData(rse_geneQuake), 
	pd[match(colnames(rse_geneQuake), pd$RunName),2:14])

rse_geneQuake$Cell_type[rse_geneQuake$Cell_type == "hybrids"] = "hybrid"
## exclude hybrid
rse_geneQuake = rse_geneQuake[,rse_geneQuake$Cell_type != "hybrid"]
##
# get expression
yExprsQuake = log2(getRPKM(rse_geneQuake, "Length")+1)

## bring in YEO iPSC data
load(file.path(MAINDIR,"yeo_singleCell/rse_gene_yeo_n214.Rdata"))

rse_geneYeo = rse_gene[,rse_gene$cell_type %in% c("iPSC", "NPC") &
	rse_gene$sample_type=="Single_Cell"]
yExprsYeo = log2(getRPKM(rse_geneYeo, "Length")+1)

## merge w/ scorecard
yExprs_Merge = cbind(yExprsYeo, yExprsQuake)	
group = c(as.character(rse_geneYeo$cell_type), 
	as.character(rse_geneQuake$Cell_type))
	
	
#####################
## analysis
##################

## drop low expression
rse_geneQuake = rse_geneQuake[rowData(rse_geneQuake)$meanExprs > 1,]

## do modeling
dge = DGEList(counts = assays(rse_geneQuake)$counts, 
	genes = rowData(rse_geneQuake))
dge = calcNormFactors(dge)

## lets just look at neurons - adult versus fetal quiescent
rse_geneQuake$isNeuron = ifelse(rse_geneQuake$Cell_type == "Neurons" , 1, 0)

vGeneNeuron = voom(dge,	model.matrix(~rse_geneQuake$isNeuron), plot=TRUE)
fitGeneNeuron = lmFit(vGeneNeuron)
eBGeneNeuron = eBayes(fitGeneNeuron)
sigGeneNeuron = topTable(eBGeneNeuron,coef=2,
	p.value = 1e-20,number=nrow(rse_geneQuake))
sigGeneNeuron$gencodeTx = NULL

## add number of 0s
zeroMat = sapply(splitit(rse_geneQuake$isNeuron), 
	function(ii) rowMeans(assays(rse_geneQuake[,ii])$counts == 0))
colnames(zeroMat) = c("prop0s_nonNeuron", "prop0s_Neuron")
sigGeneNeuron = cbind(sigGeneNeuron, zeroMat[rownames(sigGeneNeuron),])
#####
par(mar=c(8,6,2,2))
boxplot(yExprs_Merge[rownames(sigGeneNeuron)[1],] ~ group,
		las=3, main = sigGeneNeuron$Symbol[1])
		
boxplot(yExprs_Merge[rownames(sigGeneNeuron)[3],] ~ group,
	las=3, main = sigGeneNeuron$Symbol[3])
boxplot(yExprs_Merge[rownames(sigGeneNeuron)[4],] ~ group,
	las=3, main = sigGeneNeuron$Symbol[4])
boxplot(yExprs_Merge[rownames(sigGeneNeuron)[5],] ~ group,
	las=3, main = sigGeneNeuron$Symbol[5])
boxplot(yExprs_Merge[rownames(sigGeneNeuron)[8],] ~ group,
	las=3, main = sigGeneNeuron$Symbol[8])
boxplot(yExprs_Merge[rownames(sigGeneNeuron)[9],] ~ group,
	las=3, main = sigGeneNeuron$Symbol[9])