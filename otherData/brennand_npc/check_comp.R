##
library(jaffelab)
library(minfi)
library(SummarizedExperiment)
library(recount)

## read in gene
geneCounts = read.csv("geneCounts.csv", as.is=TRUE)
# pd = read.csv("info.csv",as.is=TRUE)
pd = read.csv("COS RNAseq Master Sheet - Co-variates for analysis.csv",as.is=TRUE)
geneMap = read.delim("ENSEMBLv70_gene_info.tsv",as.is=TRUE,row.names=1)

## fix name
pd$SampleID = paste0("X", gsub("-", ".", pd$Sample.Name))
pd$SampleID = gsub("FR", "F", pd$SampleID) # fix one
pd$SampleID = gsub("NR", "N", pd$SampleID) # fix two

pd$Cell.Type = ifelse(pd$Cell.Type == "NPC", "NPC", "Neuron")
pd$Cell.Type = factor(pd$Cell.Type , levels =c("NPC", "Neuron"))

# reorder
geneCounts = geneCounts[,pd$SampleID]
identical(rownames(geneCounts), rownames(geneMap)) # TRUE

## get RPKM
bg <- matrix(colSums(geneCounts), ncol = ncol(geneCounts), 
	nrow = nrow(geneCounts),byrow = TRUE)
wid <- matrix(geneMap$Length, nrow = nrow(geneCounts),
	ncol = ncol(geneCounts), byrow = FALSE)
geneRpkm = geneCounts/(wid/1000)/(bg/1e+06)

geneExprs = log2(geneRpkm+1)

###################
## load cell type
load("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/deconvolution/cell_type/singleCell_iPSC_quake_coefEsts_calibration_Zscale.rda")
yExprs_cell = scale(geneExprs[ss(rownames(coefEsts),"\\."),])
cellPropEsts = minfi:::projectCellType(yExprs_cell,coefEsts)
cellPropEsts = as.data.frame(cellPropEsts)

pd$CellDx = paste0(pd$Cell.Type, ":", pd$Dx)
pd$CellDx = factor(pd$CellDx, levels = c("NPC:CT", "NPC:SZ","Neuron:CT", "Neuron:SZ"))

gIndexes = splitit(pd$CellDx)
cellPropEsts_groupMeans = sapply(gIndexes, 
	function(ii) colMeans(cellPropEsts[ii,]))
cellPropEsts_groupSEs = sapply(gIndexes, 
	function(ii) apply(cellPropEsts[ii,],2,function(x) sd(x)/sqrt(length(x))))

pal = brewer.pal(10, "Set3")
pal[2] = "#8B795E" # darken NPCs
 
pdf("../deconvolution/cell_type/lineplots/brennand_stemcell_timecourse_cellTypeDecon.pdf",h=5,w=5)
palette(pal)
par(mar=c(10,6,2,2),cex.axis=1.2,cex.lab=1.5)
plot(cellPropEsts_groupMeans[1,1:2], type="b",xaxt="n",
	lty=1, pch=19,lwd=2, ylim = c(0,1),col=1,xlab="",
	ylab = "Class Proportion",cex=1.3,xlim = c(0.5,4.5))
axis(1,at=seq(along=gIndexes),colnames(cellPropEsts_groupMeans), las=3,cex.axis=1.5)
segments(x0=seq(along=gIndexes)[1:2], x1=seq(along=gIndexes)[1:2], col=1,lwd=2, 
	y0=cellPropEsts_groupMeans[1,1:2] - 2*cellPropEsts_groupSEs[1,1:2], 
	y1=cellPropEsts_groupMeans[1,1:2] + 2*cellPropEsts_groupSEs[1,1:2])
for(i in 2:10) {
	lines(cellPropEsts_groupMeans[i,1:2], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=seq(along=gIndexes)[1:2], x1=seq(along=gIndexes)[1:2], col=i,lwd=2, 
	y0=cellPropEsts_groupMeans[i,1:2] - 2*cellPropEsts_groupSEs[i,1:2], 
	y1=cellPropEsts_groupMeans[i,1:2] + 2*cellPropEsts_groupSEs[i,1:2])
}
for(i in 1:10) {
	lines(x=3:4, cellPropEsts_groupMeans[i,3:4], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=seq(along=gIndexes)[3:4], x1=seq(along=gIndexes)[3:4], col=i,lwd=2, 
	y0=cellPropEsts_groupMeans[i,3:4] - 2*cellPropEsts_groupSEs[i,3:4], 
	y1=cellPropEsts_groupMeans[i,3:4] + 2*cellPropEsts_groupSEs[i,3:4])
}
dev.off()

###################
## load brain stage
load("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/deconvolution/brain_stage/iPSC_BrainSpan_coefEsts_calibration_Zscore.rda")
yExprs_stage = scale(geneExprs[ss(rownames(coefEsts),"\\."),])
stagePropEsts = minfi:::projectCellType(yExprs_stage,coefEsts)
stagePropEsts = as.data.frame(stagePropEsts)

### make plots ####
mypar(2,5)
for(i in 1:ncol(cellPropEsts)) {
	boxplot(cellPropEsts[,i] ~ pd$Cell.Type,
		ylim = c(0,1), main = colnames(cellPropEsts)[i])
}
for(i in 1:ncol(cellPropEsts)) {
	boxplot(cellPropEsts[,i] ~ pd$Dx*pd$Cell.Type,
		ylim = c(0,1), main = colnames(cellPropEsts)[i])
}

lapply(cellPropEsts, function(y) {
	signif(summary(lm(y~pd$Dx*pd$Cell.Type))$coef,3)
})

boxplot(cellPropEsts$NPC ~ pd$Cell.Type)
boxplot(cellPropEsts$Fetal_replicating ~ pd$Cell.Type)
boxplot(cellPropEsts$Fetal_quiescent ~ pd$Cell.Type)
boxplot(cellPropEsts$Neurons ~ pd$Cell.Type)
boxplot(cellPropEsts$Endothelial ~ pd$Cell.Type)