###
library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(genefilter)
library(RColorBrewer)

MAINDIR = "/dcl01/lieber/ajaffe/lab/libd_stem_timecourse"


########################
### read in pheno
pheno = read.table("/dcl01/lieber/ajaffe/PublicData/Brain/Zhang_CellType/SraRunTable.txt", 
				header=TRUE, stringsAsFactors=FALSE, sep="\t")
rownames(pheno) = pheno$Run_s


########################
### full samples
load("/dcl01/lieber/ajaffe/PublicData/Brain/Zhang_CellType/hg38/rse_gene_Zhang_CellType_Human_hg38_n41.Rdata")
pd = colData(rse_gene)
stopifnot(identical(rownames(pheno), rownames(colData(rse_gene))))
pd = cbind(pd, pheno)

## No tumor samples
rse_gene = rse_gene[,-grep("tumor", pd$source_name_s)]
pd = pd[-grep("tumor", pd$source_name_s),]
## Missing info
missInd = which(pd$cell_type_s %in% c("","Neuron"))
rse_gene = rse_gene[,-missInd]
pd = pd[-missInd,]

pd$label = pd$cell_type_s
pd$label = ifelse(pd$label=="Astrocyte", paste0(pd$cell_type_s, "_", pd$source_name_s), pd$label)
pd$label = as.factor(pd$label)
pd$label = factor(pd$label, levels(pd$label)[c(1,3,2,7,4,6,5)] )

#get RPKM
yExprsTC = log2(recount::getRPKM(rse_gene,length_var="Length")+1)





########################
### deconvolution

###############################################################
############## cell type ######################################

pal = brewer.pal(10, "Set3")
pal[2] = "#8B795E" # darken NPCs

load("cell_type/singleCell_iPSC_quake_coefEsts_calibration_Zscale.rda")

yExprsTC_Z = scale(yExprsTC[rownames(coefEsts),])

# project
stemPropEsts = minfi:::projectCellType(yExprsTC_Z,coefEsts)
stemPropEstsScaled = prop.table(stemPropEsts,1)

## line plot
gIndexes_tc = splitit(pd$label)
stemPropEsts_groupMeans = sapply(gIndexes_tc, 
	function(ii) colMeans(stemPropEsts[ii,]))
stemPropEsts_groupSEs = sapply(gIndexes_tc, 
	function(ii) apply(stemPropEsts[ii,],2,function(x) sd(x)/sqrt(length(x))))

round(100*stemPropEsts_groupMeans,3)

stemPropEsts_groupMeans = prop.table(stemPropEsts_groupMeans,2)  ## add up to 100



pdf("cell_type/lineplots/zhang_astro_cellTypeDecon.pdf",h=5,w=5)
palette(pal)
par(mar=c(10,6,2,2),cex.axis=1.2,cex.lab=1.6)
plot(stemPropEsts_groupMeans[1,], type="p",xaxt="n",
	lty=1, pch=19,lwd=2, ylim = c(0,1),col=1,xlab="",
	ylab = "Class Proportion",cex=1.3)
axis(1,at=seq(along=gIndexes_tc),levels(pd$label), las=3,cex.axis=.9)
segments(x0=seq(along=gIndexes_tc), x1=seq(along=gIndexes_tc), col=1,lwd=2, 
	y0=stemPropEsts_groupMeans[1,] - 2*stemPropEsts_groupSEs[1,], 
	y1=stemPropEsts_groupMeans[1,] + 2*stemPropEsts_groupSEs[1,])
for(i in 1:10) {
	lines(stemPropEsts_groupMeans[i,1:3], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=seq(along=gIndexes_tc[1:3]), x1=seq(along=gIndexes_tc[1:3]), col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,1:3] - 2*stemPropEsts_groupSEs[i,1:3], 
	y1=stemPropEsts_groupMeans[i,1:3] + 2*stemPropEsts_groupSEs[i,1:3])
}
for(i in 1:10) {
	lines(x=4:7, stemPropEsts_groupMeans[i,4:7], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=4:7, x1=4:7, col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,4:7] - 2*stemPropEsts_groupSEs[i,4:7], 
	y1=stemPropEsts_groupMeans[i,4:7] + 2*stemPropEsts_groupSEs[i,4:7])
}
dev.off()



