###

library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(genefilter)

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

# get expression
yExprsQuake = log2(getRPKM(rse_geneQuake, "Length")+1)

## bring in scorecard data
load(file.path(MAINDIR,"scorecard/rse_gene_scorecard_hg38_n73.Rdata"))
pd = read.table(file.path(MAINDIR,"scorecard/scorecard_SraRunTable.txt"), header=TRUE, sep="\t")
stopifnot(identical(as.character(pd$Run_s), colData(rse_gene)$SAMPLE_ID))
colData(rse_gene) = colData(rse_gene)[,1:2]
colData(rse_gene)$CONDITION = pd$cell_type_s
colData(rse_gene)$Donor = pd$data_set_s
colData(rse_gene)$DX = pd$data_set_s
colData(rse_gene)$DAY = pd$data_set_s
colData(rse_gene)$LibraryBatch = pd$data_set_s
colData(rse_gene)$experiment = "scorecard"
rse_geneSC = rse_gene

rse_geneSC = rse_geneSC[,rse_geneSC$CONDITION == "iPSC"]
yExprsSC = log2(getRPKM(rse_geneSC, "Length")+1)

##########
## merge #

## merge w/ scorecard
yExprs_Merge = cbind(yExprsSC, yExprsQuake)	
group = c(as.character(rse_geneSC$CONDITION), 
	as.character(rse_geneQuake$Cell_type))
	
group = factor(group,levels =c("iPSC", "Fetal_replicating",
	"Fetal_quiescent", "OPC", "Neurons", "Astrocytes",
	"Oligodendrocytes", "Microglia", "Endothelial"))

## #split by age cat and region					
tIndexes <- splitit(group)

tstatList <- lapply(tIndexes, function(i) {
	x <- rep(0, ncol(yExprs_Merge))
	x[i] <- 1
	return(genefilter::rowttests(yExprs_Merge, factor(x)))
})
numProbes=20
probeList <- lapply(tstatList, function(x) {
	y <- x[x[, "p.value"] < 1e-15, ]
	yUp <- y[order(y[, "dm"], decreasing = TRUE), ]
	yDown <- y[order(y[, "dm"], decreasing = FALSE),]
	c(rownames(yUp)[1:numProbes], rownames(yDown)[1:numProbes])
})

## filter
trainingProbes <- unique(unlist(probeList))
mergeMarkerExprs <- yExprs_Merge[trainingProbes, ]

mergeMarkerMeanExprs <- colMeans(mergeMarkerExprs)

form <- as.formula(sprintf("y ~ %s - 1", paste(levels(group),collapse = "+")))
phenoDF <- as.data.frame(model.matrix(~group - 1))
colnames(phenoDF) <- sub("^group", "", colnames(phenoDF))

## try z-score
mergeMarkerExprsZ = scale(mergeMarkerExprs)
mergeMarkerMeanExprsZ = colMeans(mergeMarkerExprsZ)

## do calibration
coefEsts <- minfi:::validationCellType(Y = mergeMarkerExprsZ, 
	pheno = phenoDF, modelFix = form)$coefEsts

save(coefEsts, mergeMarkerMeanExprs, 
	file = "Scorecard_iPSC_quakeSingleCell_coefEsts_calibration_Zscale.rda")

#####################
### test on LIBD ####

# load( "Scorecard_iPSC_quakeSingleCell_coefEsts_calibration.rda")

load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseGene_merged_n732.rda")
rse_geneDLPFC = rse_gene[,rse_gene$Dx == "Control" & 
	rse_gene$Race %in% c("AA","CAUC")]
rse_geneDLPFC$ageGroup = cut(rse_geneDLPFC$Age, c(-1,0,1,10,20,50,100))
levels(rse_geneDLPFC$ageGroup) = c("Fetal","Infant","Child","Teens","Adult","50+")

rse_geneDLPFC = rse_geneDLPFC[,order(rse_geneDLPFC$Age)]
yExprsDLPFC = log2(getRPKM(rse_geneDLPFC,length_var="Length")+1)

# project
yExprsDLPFC_Z = scale(yExprsDLPFC[rownames(coefEsts),])
dlpfcPropEsts = minfi:::projectCellType(yExprsDLPFC_Z,coefEsts)
dlpfcPropEstsScaled = prop.table(dlpfcPropEsts,1)
pal = brewer.pal(9,"Set1")

## dlpfc
pdf("Quake/brainSeq_compEsts.pdf", w=14,h=6)
palette(pal)
par(mar = c(6,6,2,2), cex.axis=1.5,cex.lab=2)
ooDlpfc = order(rse_geneDLPFC$Age)
bp=barplot(t(dlpfcPropEstsScaled[ooDlpfc,]), col = pal,
	ylim = c(0,1.4),xaxt = "n", yaxt= "n",ylab="Class Proportion")

g = split(bp, rse_geneDLPFC$ageGroup[ooDlpfc])
axis(1,at=sapply(g,mean), font=2, labels = names(g),las=3)
abline(v = sapply(g,min)[-1]-0.5,lwd=4,col="black")
legend("top", colnames(dlpfcPropEstsScaled), pch = 15, 
	col = 1:9,bg="white", nc = 4,cex=1.5,pt.cex=2)
axis(2, at=seq(0,1,by=0.25))
dev.off()

## diagnose
plot(c(mergeMarkerMeanExprsZ, colMeans(yExprsDLPFC_Z)),
		pch = 21, bg = c(as.numeric(group), 
			rep("black", ncol(yExprsDLPFC))),
	ylab = "meanExprs")

##################
# libd stem ######
## libd time course
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n146.rda"))
getRPKM = recount::getRPKM

## don't use technical replicates
bioInd = which(colData(rse_gene)$Class == "Naked genomes")
rse_gene = rse_gene[,bioInd]  # n=128
colData(rse_gene) = colData(rse_gene)[,c(24,26,23,15,17,20,6)]
colData(rse_gene)$DAY = paste0("day_",colData(rse_gene)$DAY)
colData(rse_gene)$experiment = "lieber"
rse_geneTC = rse_gene

## clean CONDITION
rse_geneTC$COND = rse_geneTC$CONDITION
rse_geneTC$COND[grep("ACC_DORSAL",rse_geneTC$COND)] = "ACC_DORSAL"
rse_geneTC$COND[grep("PLUS",rse_geneTC$COND)] = "NEURONS+ASTROS"
rse_geneTC$COND = as.factor(rse_geneTC$COND)
rse_geneTC$COND = factor(rse_geneTC$COND,levels(rse_geneTC$COND)[c(5,1,4,6,2,3)])
rse_geneTC$DAYNUM = as.numeric(ss(rse_geneTC$DAY, "_",2))

#get RPKM
yExprsTC = log2(getRPKM(rse_geneTC,length_var="Length")+1)
yExprsTC_Z = scale(yExprsTC[rownames(coefEsts),])

# project
stemPropEsts = minfi:::projectCellType(yExprsTC_Z,coefEsts)
stemPropEstsScaled = prop.table(stemPropEsts,1)

pdf("Quake/LIBDstem_compEsts.pdf",w=14,h=8)
palette(pal)
par(mar = c(14,6,2,2), cex.axis=1.5,cex.lab=2)
ooStem = order(rse_geneTC$COND, rse_geneTC$DAYNUM)
bp = barplot(t(stemPropEstsScaled[ooStem,]), col = pal,
	ylim = c(0,1.4),xaxt = "n", yaxt= "n",ylab="Class Proportion")
g = split(bp, rse_geneTC$COND[ooStem])

axis(1,at=sapply(g,mean), font=2, labels = names(g),las=3)
abline(v = sapply(g,min)[-1]-0.5,lwd=2,col="black")
legend("top", colnames(stemPropEstsScaled), pch = 15, 
	col = 1:9,bg="white", nc = 4,cex=1.5,pt.cex=2)
axis(2, at=seq(0,1,by=0.25))
dev.off()
