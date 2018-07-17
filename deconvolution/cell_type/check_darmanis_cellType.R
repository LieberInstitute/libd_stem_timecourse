###

library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(genefilter)
library(RColorBrewer)

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

# get expression
yExprsQuake = log2(getRPKM(rse_geneQuake, "Length")+1)

## bring in YEO iPSC data
load(file.path(MAINDIR,"yeo_singleCell/rse_gene_yeo_n214.Rdata"))

rse_geneYeo = rse_gene[,rse_gene$cell_type %in% c("iPSC", "NPC") &
	rse_gene$sample_type=="Single_Cell"]
yExprsYeo = log2(getRPKM(rse_geneYeo, "Length")+1)

##########
## merge #

## merge w/ scorecard
yExprs_Merge = cbind(yExprsYeo, yExprsQuake)	
group = c(as.character(rse_geneYeo$cell_type), 
	as.character(rse_geneQuake$Cell_type))
	
group = factor(group,levels =c("iPSC", "NPC","Fetal_replicating",
	"Fetal_quiescent", "OPC", "Neurons", "Astrocytes",
	"Oligodendrocytes", "Microglia", "Endothelial"))

## #split by age cat and region					
tIndexes <- splitit(group)

tstatList <- lapply(tIndexes, function(i) {
	x <- rep(0, ncol(yExprs_Merge))
	x[i] <- 1
	return(genefilter::rowttests(yExprs_Merge, factor(x)))
})

numProbes=25
probeList <- lapply(tstatList, function(x) {
	y <- x[which(x[, "p.value"] < 1e-15), ]
	yUp <- y[order(y[, "dm"], decreasing = FALSE),] # signs are swapped
	rownames(yUp)[1:numProbes]
})

## filter
trainingProbes <- unique(unlist(probeList))
trainingProbes = trainingProbes[!is.na(trainingProbes)]

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

save(coefEsts, mergeMarkerMeanExprsZ, 
	file = "singleCell_iPSC_quake_coefEsts_calibration_Zscale.rda")
write.csv(coefEsts, file = "singleCell_iPSC_quake_coefEsts_calibration_Zscale.csv")

## add symbol
tabOut = as.data.frame(coefEsts)
tabOut$MarkerClass =  colnames(coefEsts)[apply(coefEsts, 1, which.max)]
tabOut$Symbol = rowData(rse_geneQuake[rownames(tabOut),])$Symbol
tabOut$GencodeID = rownames(tabOut)
tabOut = tabOut[,c(13,12,11,1:10)]
write.csv(tabOut, file="annotated_singleCell_iPSC_quake_coefEsts_calibration_Zscale.csv",
	row.names=FALSE)

## make plot
library(lattice)
theSeq = seq(-3.5,3.5,by=0.2)
mat =  coefEsts
colnames(mat) = c("iPSC", "NPC", "Fetal:Repl", "Fetal:Quies", "Adult:OPC",
	"Adult:Neuron", "Adult:Astro", "Adult:Oligo", "Adult:Microglia", "Adult:Endothelial")
rownames(mat) = rowData(rse_geneQuake[rownames(coefEsts),])$Symbol
	my.col <- colorRampPalette(c("blue","white","red"))(length(theSeq))
pdf("singleCellGroup_exprsMatZ.pdf",w=36)
print(levelplot(mat, aspect = "fill", at = theSeq,pretty=TRUE,
	panel = panel.levelplot.raster, col.regions = my.col,
		scales=list(x=list(rot=90,cex=1.15), y=list(cex=1.4)),
		ylab = "Cell Type", xlab = ""))
dev.off()

## drop 4 classes
mainIndex = which(tabOut$MarkerClass %in% c("Endothelial",
	"Fetal_quiescent", "Fetal_replicating", "iPSC", "Neurons", "NPC"))

pdf("singleCellGroup_exprsMatZ_main.pdf",w=36)
print(levelplot(mat[mainIndex,], aspect = "fill", at = theSeq,pretty=TRUE,
	panel = panel.levelplot.raster, col.regions = my.col,
		scales=list(x=list(rot=90,cex=2), y=list(cex=1.4)),
		ylab = "Cell Type", xlab = ""))
dev.off()

## actual raw expression
print(levelplot(yExprs_Merge[rownames(coefEsts),],
	aspect = "fill", pretty=TRUE,
	panel = panel.levelplot.raster, 
		scales=list(x=list(rot=90)),ylab = "Cell Type", xlab = ""))
	
library(pheatmap)
eMat =yExprs_Merge[rownames(coefEsts),]
colnames(eMat) = colnames(mat)[match(group, colnames(coefEsts))]
rownames(eMat) = rowData(rse_geneQuake[rownames(coefEsts),])$Symbol
pdf("singleCell_exprsMat.pdf", h=30,w=70)
print(pheatmap(eMat,
	color = colorRampPalette(brewer.pal(n = 7, name = "PuRd"))(100)))
dev.off()
		
#####################
### test on LIBD ####

load( "singleCell_iPSC_quake_coefEsts_calibration_Zscale.rda")

load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseGene_merged_n732.rda")
rse_geneDLPFC = rse_gene[,rse_gene$Dx == "Control" & 
	rse_gene$Race %in% c("AA","CAUC")]
rse_geneDLPFC$ageGroup = cut(rse_geneDLPFC$Age, c(-1,-0.45, 0,1,10,20,100))
levels(rse_geneDLPFC$ageGroup) = c("EarlyFetal","MidFetal",
	"Infant","Child","Teens","Adult")

rse_geneDLPFC = rse_geneDLPFC[,order(rse_geneDLPFC$Age)]
yExprsDLPFC = log2(getRPKM(rse_geneDLPFC,length_var="Length")+1)

# project
yExprsDLPFC_Z = scale(yExprsDLPFC[rownames(coefEsts),])
dlpfcPropEsts = minfi:::projectCellType(yExprsDLPFC_Z,coefEsts)
dlpfcPropEstsScaled = prop.table(dlpfcPropEsts,1)

##############
## make line plots

# summarize
gIndexes_dlpfc = splitit(rse_geneDLPFC$ageGroup)
dlpfcPropEsts_groupMeans = sapply(gIndexes_dlpfc, 
	function(ii) colMeans(dlpfcPropEsts[ii,]))
dlpfcPropEsts_groupSEs = sapply(gIndexes_dlpfc, 
	function(ii) apply(dlpfcPropEsts[ii,],2,function(x) sd(x)/sqrt(length(x))))

## make legend pdf
pdf("lineplots/legend_cellTypeDecon.pdf")
plot(0,0, axes=FALSE,type="n")
legend("top", rownames(dlpfcPropEsts_groupMeans), 
	ncol = 1, pch = 15, col = 1:10,pt.cex=1.5,cex=1.3)
dev.off()


## make line plot
pdf("lineplots/brainseq_phase1_dlpfc_cellTypeDecon.pdf",h=5,w=5)
palette(pal)
par(mar=c(10,6,2,2),cex.axis=1.2,cex.lab=1.5)
plot(dlpfcPropEsts_groupMeans[1,], type="b",xaxt="n",
	lty=1, pch=19,lwd=2, ylim = c(0,1),col=1,xlab="",
	ylab = "Class Proportion",cex=1.3)
axis(1,at=seq(along=gIndexes_dlpfc),colnames(dlpfcPropEsts_groupMeans), 
	las=3,cex.axis=1.5)
segments(x0=seq(along=gIndexes_dlpfc), x1=seq(along=gIndexes_dlpfc), col=1,lwd=2, 
	y0=dlpfcPropEsts_groupMeans[1,] - 2*dlpfcPropEsts_groupSEs[1,], 
	y1=dlpfcPropEsts_groupMeans[1,] + 2*dlpfcPropEsts_groupSEs[1,])

for(i in 2:10) {
	lines(dlpfcPropEsts_groupMeans[i,], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=seq(along=gIndexes_dlpfc), x1=seq(along=gIndexes_dlpfc), col=i,lwd=2, 
	y0=dlpfcPropEsts_groupMeans[i,] - 2*dlpfcPropEsts_groupSEs[i,], 
	y1=dlpfcPropEsts_groupMeans[i,] + 2*dlpfcPropEsts_groupSEs[i,])
}
abline(v=2.5,lty=2,col="grey")
dev.off()	

## dlpfc
pdf("barplots/brainSeq_compEsts_viaSingleCell.pdf", w=14,h=6)
palette(pal)
par(mar = c(6,6,2,2), cex.axis=1.5,cex.lab=2)
ooDlpfc = order(rse_geneDLPFC$Age)
bp=barplot(t(dlpfcPropEstsScaled[ooDlpfc,]), col = pal,
	ylim = c(0,1.4),xaxt = "n", yaxt= "n",ylab="Class Proportion")

g = split(bp, rse_geneDLPFC$ageGroup[ooDlpfc])
axis(1,at=sapply(g,mean), font=2, labels = names(g),las=3)
abline(v = sapply(g,min)[-1]-0.5,lwd=4,col="black")
legend("top", colnames(dlpfcPropEstsScaled), pch = 15, 
	col = 1:10,bg="white", nc = 5,cex=1.5,pt.cex=2)
axis(2, at=seq(0,1,by=0.25))
dev.off()

## diagnose
plot(c(mergeMarkerMeanExprsZ, colMeans(yExprsDLPFC_Z)),
		pch = 21, bg = c(as.numeric(group), 
			rep("black", ncol(yExprsDLPFC))),
	ylab = "meanExprs")
save(dlpfcPropEsts, file= "brainSeq_compEsts_viaSingleCell.rda")

## regression
signif(100*t(apply(dlpfcPropEsts,2,function(y) {
	summary(lm(y ~ as.numeric(rse_geneDLPFC$ageGroup)))$coef[2,]
})),3)

##################
# libd stem ######
## libd time course

load("singleCell_iPSC_quake_coefEsts_calibration_Zscale.rda")

###########
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n146.rda"))
getRPKM = recount::getRPKM

## don't use technical replicates
bioInd = which(colData(rse_gene)$Class == "Naked genomes")
rse_gene = rse_gene[,bioInd]  # n=128
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

###
modAdd = model.matrix(~ as.numeric(COND) + LINE + totalAssignedGene, 
	data=colData(rse_geneTC))
signif(t(apply(stemPropEsts, 2, function(y) summary(lm(y ~ modAdd - 1))$coef[2,])),3)       
	
## line plot
gIndexes_tc = splitit(rse_geneTC$COND)
stemPropEsts_groupMeans = sapply(gIndexes_tc, 
	function(ii) colMeans(stemPropEsts[ii,]))
stemPropEsts_groupSEs = sapply(gIndexes_tc, 
	function(ii) apply(stemPropEsts[ii,],2,function(x) sd(x)/sqrt(length(x))))

round(100*stemPropEsts_groupMeans,3)

pdf("lineplots/libd_stemcell_timecourse_cellTypeDecon.pdf",h=5,w=5)
palette(pal)
par(mar=c(10,6,2,2),cex.axis=1.2,cex.lab=1.5)
plot(stemPropEsts_groupMeans[1,], type="b",xaxt="n",
	lty=1, pch=19,lwd=2, ylim = c(0,1),col=1,xlab="",
	ylab = "Class Proportion",cex=1.3)
axis(1,at=seq(along=gIndexes_tc),c("RENEW","ACCEL\nDORSAL", "NPC", "ROSETTE", 
	"NEURONS","NEURONS\n+ASTROS"), las=3,cex.axis=1.5)
segments(x0=seq(along=gIndexes_tc), x1=seq(along=gIndexes_tc), col=1,lwd=2, 
	y0=stemPropEsts_groupMeans[1,] - 2*stemPropEsts_groupSEs[1,], 
	y1=stemPropEsts_groupMeans[1,] + 2*stemPropEsts_groupSEs[1,])
for(i in 2:10) {
	lines(stemPropEsts_groupMeans[i,], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=seq(along=gIndexes_tc), x1=seq(along=gIndexes_tc), col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,] - 2*stemPropEsts_groupSEs[i,], 
	y1=stemPropEsts_groupMeans[i,] + 2*stemPropEsts_groupSEs[i,])
}
dev.off()

## barplot
neuronIndex =grep("^NEURON", rse_geneTC$CONDITION )
pdf("barplots/LIBDstem_compEsts_viaSingleCell_justNeurons.pdf",w=4,h=5)
palette(pal)
par(mar=c(10,6,2,2),cex.axis=1.2,cex.lab=1.5)
bp = barplot(t(stemPropEstsScaled[neuronIndex,]), col = pal,
	ylim = c(0,1),xaxt = "n", yaxt= "n",ylab="Class Proportion")
g = split(bp, rse_geneTC$COND[neuronIndex])

axis(1,at=sapply(g,mean)[5:6], font=2, labels = c("",""),las=3)
abline(v = sapply(g,min)[6]-0.5,lwd=2,col="black")
axis(2, at=seq(0,1,by=0.25))
dev.off()


pdf("barplots/LIBDstem_compEsts_viaSingleCell.pdf",w=14,h=8)
palette(pal)
par(mar = c(14,6,2,2), cex.axis=1.5,cex.lab=2)
ooStem = order(rse_geneTC$COND, rse_geneTC$DAYNUM)
bp = barplot(t(stemPropEstsScaled[ooStem,]), col = pal,
	ylim = c(0,1.3),xaxt = "n", yaxt= "n",ylab="Class Proportion")
g = split(bp, rse_geneTC$COND[ooStem])

axis(1,at=sapply(g,mean), font=2, labels = names(g),las=3)
abline(v = sapply(g,min)[-1]-0.5,lwd=2,col="black")
legend("top", colnames(stemPropEstsScaled), pch = 15, 
	col = 1:10,bg="white", nc = 5,cex=1.5,pt.cex=2)
axis(2, at=seq(0,1,by=0.25))
dev.off()

pdf("LIBDstem_compEsts_viaSingleCell_wide.pdf",w=25,h=8)
palette(pal)
par(mar = c(14,6,2,2), cex.axis=1.5,cex.lab=2)
ooStem = order(rse_geneTC$COND, rse_geneTC$DAYNUM)
bp = barplot(t(stemPropEstsScaled[ooStem,]), col = pal,
	ylim = c(0,1.3),xaxt = "n", yaxt= "n",ylab="Class Proportion")
g = split(bp, rse_geneTC$COND[ooStem])

axis(1,at=sapply(g,mean), font=2, labels = names(g),las=3)
abline(v = sapply(g,min)[-1]-0.5,lwd=2,col="black")
legend("top", colnames(stemPropEstsScaled), pch = 15, 
	col = 1:10,bg="white", nc = 5,cex=1.5,pt.cex=2)
axis(2, at=seq(0,1,by=0.25))
dev.off()

### heatmap of signature ###

## make plot
library(lattice)
theSeq = seq(-3.5,3.5,by=0.2)
eMat =yExprsTC_Z
colnames(eMat) = rse_geneTC$COND
rownames(eMat) = rowData(rse_geneTC[rownames(coefEsts),])$Symbol
my.col <- colorRampPalette(c("blue","white","red"))(length(theSeq))

pdf("humanNeuron_exprsMatZ_main.pdf",w=36,h=10)
print(levelplot(eMat[,oo], aspect = "fill", at = theSeq,pretty=TRUE,
	panel = panel.levelplot.raster, col.regions = my.col,
		scales=list(x=list(rot=90,cex=1.15)),	
		ylab = "", xlab = ""))
dev.off()

pdf("humanNeuron_exprsMatZ_onlyMain.pdf",w=36,h=10)
print(levelplot(eMat[mainIndex,ooStem], aspect = "fill", at = theSeq,pretty=TRUE,
	panel = panel.levelplot.raster, col.regions = my.col,
		scales=list(x=list(rot=90,cex=1.85)),	
		ylab = "", xlab = ""))
dev.off()

########################
## add symbol to blank ones (per ensembl.org)
# ensemblID
# ENSG00000254277	AC009446.1
# ENSG00000249152	LNCPRESS2
# ENSG00000280707	AL353747.4
# ENSG00000254339	AC064802.1
# ENSG00000253507	AC104257.1
# ENSG00000280511	AL591030.1
# ENSG00000256637	LINC01965
# ENSG00000278543	AC010332.2
# ENSG00000201774	RF00019
# ENSG00000253638	AF186190.1
# ENSG00000280234	AC124303.2
syms = c("AC009446.1","LNCPRESS2","AL353747.4","AC064802.1","AC104257.1",
			"AL591030.1","LINC01965","AC010332.2","RF00019","AF186190.1","AC124303.2")

ooStem = order(rse_geneTC$COND, rse_geneTC$DAYNUM)
eMat_sub = eMat[mainIndex,ooStem]
eMap_sub = rowData(rse_geneTC[rownames(coefEsts),])
eMap_sub = eMap_sub[mainIndex,]
eMap_sub$Symbol[which(eMap_sub$Symbol=="")] = syms
rownames(eMat_sub) = eMap_sub$Symbol

pdf("humanNeuron_exprsMatZ_onlyMain_eeb.pdf",w=36,h=10)
print(levelplot(eMat_sub, aspect = "fill", at = theSeq,pretty=TRUE,
	panel = panel.levelplot.raster, col.regions = my.col,
		scales=list(x=list(rot=90,cex=1.85)),	
		ylab = "", xlab = ""))
dev.off()


## make plot of weights
library(lattice)
theSeq2 = seq(0,1,by=0.025)
wts = stemPropEstsScaled
colnames(wts) = colnames(mat)
rownames(wts) = rse_geneTC$COND
my.col2 <- colorRampPalette(c("white","black"))(length(theSeq2))
oo2 = order(colMeans(stemPropEstsScaled), decreasing=TRUE)
pdf("humanNeuron_cellCompWeights.pdf")
print(levelplot(wts[ooStem,oo2], aspect = "fill",
	at = theSeq2,pretty=TRUE,
	panel = panel.levelplot.raster, col.regions = my.col2,
		scales=list(x=list(rot=90,cex=0.5)),	
		ylab = "", xlab = ""))
dev.off()


library(pheatmap)

pdf("humanNeuron_exprsMat.pdf", h=30,w=70)
print(pheatmap(eMat,
	color = colorRampPalette(brewer.pal(n = 7, name = "PuRd"))(100)))
dev.off()

#### examples of genes
reference = scale(yExprs_Merge[rownames(coefEsts),])
sym = rowData(rse_geneTC)$Symbol

g = c("SCG2", "NEUROD6")
palette(pal)

boxplot(yExprs_Merge[match("NEUROD6", rowData(rse_geneTC)$Symbol),] ~ group)
boxplot(reference[match("NEUROD6", sym),] ~ group)
boxplot(yExprs_Merge[match("SCG2", rowData(rse_geneTC)$Symbol),] ~ group)
boxplot(reference[match("SCG2", sym),] ~ group)

boxplot(yExprs_Merge[match("SCG2", rowData(rse_geneTC)$Symbol),] ~ group)

plot(reference[match("NEUROD6", sym),], reference[match("SCG2", sym),],
	pch = 21, bg = group)

##############################
## ALL THE DATA!!! ###########
##############################

######################
## load cortecon data ## hg38
load("/dcs01/ajaffe/LIBD_StemCell/Data/processed_rawCounts_Cortecon.rda") ## for Day, other pd
pdCort2 = pd[,1:8]
load("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/cortecon/rpkmCounts_cortecon_hg38_n24.rda")
load("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/cortecon/rawCounts_cortecon_hg38_n24.rda")
pdCort = metrics
stopifnot(identical(rownames(pdCort),pdCort2$SRR))
pdCort = cbind(pdCort,pdCort2)
geneCountsCort = geneCounts
geneRpkmCort = geneRpkm
geneMapCort = geneMap
rm(exonCounts, exonRpkm, exonMap, jCounts, jRpkm, jMap, txNumReads, txTpm, txMap, pdCort2)

yExprsCort = log2(geneRpkmCort+1)
yExprsCort_Scaled = scale(yExprsCort[rownames(coefEsts),])

# project
cortPropEsts = minfi:::projectCellType(yExprsCort_Scaled,coefEsts)
cortPropEstsScaled = prop.table(cortPropEsts,1)

# line plot
gIndexes_cort = splitit(pdCort$Day)
cortPropEsts_groupMeans = sapply(gIndexes_cort, 
	function(ii) colMeans(cortPropEsts[ii,]))
cortPropEsts_groupSEs = sapply(gIndexes_cort, 
	function(ii) apply(cortPropEsts[ii,],2,function(x) sd(x)/sqrt(length(x))))

pdf("lineplots/cortecon_cellTypeDecon.pdf",h=5,w=5)
palette(pal)
par(mar=c(10,6,2,2),cex.axis=1.2,cex.lab=1.3)
plot(cortPropEsts_groupMeans[1,], type="b",xaxt="n",
	lty=1, pch=19,lwd=2, ylim = c(0,1),col=1,xlab="",
	ylab = "Class Proportion",cex=1.3)
axis(1,at=seq(along=gIndexes_cort),gsub(" ", "", colnames(cortPropEsts_groupMeans)), 
	las=3,cex.axis=1.5)
segments(x0=seq(along=gIndexes_cort), x1=seq(along=gIndexes_cort), col=1,lwd=2, 
	y0=cortPropEsts_groupMeans[1,] - 2*cortPropEsts_groupSEs[1,], 
	y1=cortPropEsts_groupMeans[1,] + 2*cortPropEsts_groupSEs[1,])

for(i in 2:10) {
	lines(cortPropEsts_groupMeans[i,], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=seq(along=gIndexes_cort), x1=seq(along=gIndexes_cort), col=i,lwd=2, 
	y0=cortPropEsts_groupMeans[i,] - 2*cortPropEsts_groupSEs[i,], 
	y1=cortPropEsts_groupMeans[i,] + 2*cortPropEsts_groupSEs[i,])
}
dev.off()	

## human data from span, homogenate
load("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/deconvolution/brain_stage/brainspan_gene_rse.rda")
yExprsSpan = log2(getRPKM(rse_geneSpan,length_var="Length")+1)
yExprsSpan_Scaled = scale(yExprsSpan[rownames(coefEsts),])

# project
spanPropEsts = minfi:::projectCellType(yExprsSpan_Scaled,coefEsts)
spanPropEstsScaled = prop.table(spanPropEsts,1)

# line plot
gIndexes_span = splitit(rse_geneSpan$AgeCat)
spanPropEsts_groupMeans = sapply(gIndexes_span, 
	function(ii) colMeans(spanPropEsts[ii,]))
spanPropEsts_groupSEs = sapply(gIndexes_span, 
	function(ii) apply(spanPropEsts[ii,],2,function(x) sd(x)/sqrt(length(x))))

pdf("lineplots/brainspan_homogenate_cellTypeDecon.pdf",h=5,w=5)
palette(pal)
par(mar=c(10,6,2,2),cex.axis=1.2,cex.lab=1.3)
plot(spanPropEsts_groupMeans[1,], type="b",xaxt="n",
	lty=1, pch=19,lwd=2, ylim = c(0,1),col=1,xlab="",
	ylab = "Class Proportion",cex=1.3)
axis(1,at=seq(along=gIndexes_span),gsub(" ", "", colnames(spanPropEsts_groupMeans)), 
	las=3,cex.axis=1.5)
segments(x0=seq(along=gIndexes_span), x1=seq(along=gIndexes_span), col=1,lwd=2, 
	y0=spanPropEsts_groupMeans[1,] - 2*spanPropEsts_groupSEs[1,], 
	y1=spanPropEsts_groupMeans[1,] + 2*spanPropEsts_groupSEs[1,])

for(i in 2:10) {
	lines(spanPropEsts_groupMeans[i,], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=seq(along=gIndexes_span), x1=seq(along=gIndexes_span), col=i,lwd=2, 
	y0=spanPropEsts_groupMeans[i,] - 2*spanPropEsts_groupSEs[i,], 
	y1=spanPropEsts_groupMeans[i,] + 2*spanPropEsts_groupSEs[i,])
}
dev.off()	

pdf("barplots/BrainSpan_homogenate_compEsts_viaSingleCell.pdf",w=20,h=8)
palette(pal)
par(mar = c(14,6,2,2), cex.axis=1.5,cex.lab=2)
ooSpan= order(rse_geneSpan$Age, rse_geneSpan$Regioncode)
bp = barplot(t(spanPropEstsScaled[ooSpan,]), col = pal,
	ylim = c(0,1.4),xaxt = "n", yaxt= "n",ylab="Class Proportion")
g = split(bp, rse_geneSpan$AgeCat[ooSpan])

axis(1,at=sapply(g,mean), font=2, labels = names(g),las=3)
abline(v = sapply(g,min)[-1]-0.5,lwd=2,col="black")
legend("top", colnames(spanPropEstsScaled), pch = 15, 
	col = 1:10,bg="white", nc = 5,cex=1.5,pt.cex=2)
axis(2, at=seq(0,1,by=0.25))
dev.off()

#################
## load brainspan single cell
load("/dcl01/lieber/ajaffe/psychENCODE_Data/Sestan_SingleCell/rse_gene_brainSpan_singleCell_hg38_n932.Rdata")
rse_geneSpanSingle = rse_gene
pdSpanSingle = read.csv("/dcl01/lieber/ajaffe/psychENCODE_Data/Sestan_SingleCell/Yale_P50MH106934_SingleCellRNAseq_Metadata_RNAseq_August2016Release.csv",
	as.is=TRUE)
phenoSpanSingle = read.csv("/dcl01/lieber/ajaffe/psychENCODE_Data/Sestan_SingleCell/Yale_P50MH106934_SingleCellRNAseq_Metadata_Clinical_August2016Release.csv",
	as.is=TRUE)
phenoSpanSingle$PCW = as.numeric(ss(phenoSpanSingle$AgeDeath, "P"))
pdSpanSingle$PCW = phenoSpanSingle$PCW[match(pdSpanSingle$Individual_ID, 
	phenoSpanSingle$Individual_ID)]
rownames(pdSpanSingle) = paste0(pdSpanSingle$Individual_ID, "_",
	ss(pdSpanSingle$Sample_ID, "\\.",2), "_",
	ss(pdSpanSingle$Sample_ID, "\\.",3))
pdSpanSingle = pdSpanSingle[colnames(rse_geneSpanSingle),]
rse_geneSpanSingle$PCW = pdSpanSingle$PCW

## expression filter
rse_geneSpanSingle = rse_geneSpanSingle[,rse_geneSpanSingle$numReads > 1e5 ]
yExprs_SpanSingle = log2(getRPKM(rse_geneSpanSingle,"Length")+1)
yExprs_SpanSingle_Scaled = scale(yExprs_SpanSingle[rownames(coefEsts),])

spanSingle_PropEsts= minfi:::projectCellType(yExprs_SpanSingle[rownames(coefEsts),],coefEsts)
spanSingle_PropEstsScaled = prop.table(spanSingle_PropEsts,1)

## line plot
gIndexes_spanSingle = splitit(rse_geneSpanSingle$PCW)
spanSingle_PropEsts_groupMeans = sapply(gIndexes_spanSingle, 
	function(ii) colMeans(spanSingle_PropEstsScaled[ii,]))
spanSingle_PropEsts_groupSEs = sapply(gIndexes_spanSingle, 
	function(ii) apply(spanSingle_PropEstsScaled[ii,],2,function(x) sd(x)/sqrt(length(x))))

pdf("lineplots/brainspan_singleCell_cellTypeDecon.pdf",h=5,w=5)
palette(pal)
par(mar=c(10,6,2,2),cex.axis=1.2,cex.lab=1.3)
plot(spanSingle_PropEsts_groupMeans[1,], type="b",xaxt="n",
	lty=1, pch=19,lwd=2, ylim = c(0,1),col=1,xlab="",
	ylab = "Class Proportion",cex=1.3)
axis(1,at=seq(along=gIndexes_spanSingle),
	paste0(colnames(spanSingle_PropEsts_groupMeans),"_PCW"), las=3,cex.axis=1.5)
segments(x0=seq(along=gIndexes_spanSingle), 
	x1=seq(along=gIndexes_spanSingle),col=1,lwd=2, 
	y0=spanSingle_PropEsts_groupMeans[1,] - 2*spanSingle_PropEsts_groupSEs[1,], 
	y1=spanSingle_PropEsts_groupMeans[1,] + 2*spanSingle_PropEsts_groupSEs[1,])

for(i in 2:10) {
	lines(spanSingle_PropEsts_groupMeans[i,], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=seq(along=gIndexes_spanSingle), 
		x1=seq(along=gIndexes_spanSingle), col=i,lwd=2, 
	y0=spanSingle_PropEsts_groupMeans[i,] - 2*spanSingle_PropEsts_groupSEs[i,], 
	y1=spanSingle_PropEsts_groupMeans[i,] + 2*spanSingle_PropEsts_groupSEs[i,])
}
dev.off()	

## brainspan single cell
pdf("barplots/brainSpanSingle_compEsts_viaSingleCell.pdf", w=30,h=6)
palette(pal)
par(mar = c(6,6,2,2), cex.axis=1.5,cex.lab=2)
ooSpanSingle = order(rse_geneSpanSingle$PCW)
bp=barplot(t(spanSingle_PropEstsScaled[ooSpanSingle,]), col = pal,
	ylim = c(0,1.3),xaxt = "n", yaxt= "n",ylab="Class Proportion")

g = split(bp,rse_geneSpanSingle$PCW[ooSpanSingle])
axis(1,at=sapply(g,mean), font=2, labels = names(g),las=3)
abline(v = sapply(g,min)[-1]-0.5,lwd=4,col="black")
legend("top", colnames(spanSingle_PropEstsScaled), pch = 15, 
	col = 1:10,bg="white", nc = 5,cex=1.5,pt.cex=2)
axis(2, at=seq(0,1,by=0.25))
dev.off()

## diagnosis
maxClass = apply(spanSingle_PropEstsScaled, 1, max)
boxplot(maxClass ~ rse_geneSpanSingle$PCW)
summary(lm(maxClass ~ rse_geneSpanSingle$PCW))

################
## close iPSC ##

# load single cell data
load("/dcl01/ajaffe/data/lab/singleCell/close_interneurons/rse_gene_Close_interneurons_human_n1773.Rdata")

## add pheno
pheno = read.delim("/dcl01/ajaffe/data/lab/singleCell/close_interneurons/close_SraRunTable.txt",as.is=TRUE)
pheno = pheno[match(colnames(rse_gene), pheno$Run_s),]
colData(rse_gene) = cbind(colData(rse_gene), pheno)
rse_gene$DIV = as.numeric(gsub("D","", rse_gene$days_in_culture_s))
rse_gene_close = rse_gene
rse_gene_close$experiment = "close"
rse_gene_close$Type = ifelse(rse_gene_close$source_name_s == "cultured embryonic stem cells",
	"Bulk", "SingleCell")
# rse_gene_close$SampleLabel = paste0(rse_gene_close$cre_line_s, "_", rse_gene_close$DIV)
rse_gene_close$SampleLabel = paste0(rse_gene_close$cre_line_s, ": day ",rse_gene_close$DIV)
rse_gene_close$SampleLabel = as.factor(rse_gene_close$SampleLabel)
rse_gene_close$SampleLabel = factor(as.character(rse_gene_close$SampleLabel), 
	levels = levels(rse_gene_close$SampleLabel)[c(3,7,4,8,1,5,2,6)])

yExprs_close = log2(getRPKM(rse_gene_close,"Length")+1)
yExprs_close_Scaled = scale(yExprs_close[rownames(coefEsts),])

close_PropEsts= minfi:::projectCellType(yExprs_close_Scaled, coefEsts)
close_PropEsts_Scaled = prop.table(close_PropEsts,1)

## line plot
gIndexes_close = splitit(rse_gene_close$SampleLabel)
close_PropEsts_groupMeans = sapply(gIndexes_close, 
	function(ii) colMeans(close_PropEsts_Scaled[ii,]))
close_PropEsts_groupSEs = sapply(gIndexes_close, 
	function(ii) apply(close_PropEsts_Scaled[ii,],2,function(x) sd(x)/sqrt(length(x))))

pdf("lineplots/close_singleCell_cellTypeDecon.pdf",h=5,w=5)
palette(pal)
par(mar=c(10,6,2,2),cex.axis=1.2,cex.lab=1.3)

## DCX-
plot(close_PropEsts_groupMeans[1,seq(1,8,2)], type="b",xaxt="n",
	lty=1, pch=19,lwd=2, ylim = c(0,1),col=1,xlab="",
	ylab = "Class Proportion",cex=1.3,xlim=c(0.5,8.5))
axis(1,at=1:8, gsub(": day ", ":d", 
	colnames(close_PropEsts_groupMeans)[c(seq(1,8,2),seq(2,8,2))]), las=3,cex.axis=1.5)
segments(x0=1:4, x1=1:4,col=1,lwd=2, 
	y0=close_PropEsts_groupMeans[1,seq(1,8,2)] - 2*close_PropEsts_groupSEs[1,seq(1,8,2)], 
	y1=close_PropEsts_groupMeans[1,seq(1,8,2)] + 2*close_PropEsts_groupSEs[1,seq(1,8,2)])
for(i in 2:10) {
	lines(close_PropEsts_groupMeans[i,seq(1,8,2)], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=1:4,x1=1:4,col=i,lwd=2, 
	y0=close_PropEsts_groupMeans[i,seq(1,8,2)] - 2*close_PropEsts_groupSEs[i,seq(1,8,2)], 
	y1=close_PropEsts_groupMeans[i,seq(1,8,2)] + 2*close_PropEsts_groupSEs[i,seq(1,8,2)])
}
## DCX+
for(i in 1:10) {
	lines(x = 5:8, close_PropEsts_groupMeans[i,seq(2,8,2)], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=5:8,x1=5:8,col=i,lwd=2, 
	y0=close_PropEsts_groupMeans[i,seq(2,8,2)] - 2*close_PropEsts_groupSEs[i,seq(2,8,2)], 
	y1=close_PropEsts_groupMeans[i,seq(2,8,2)] + 2*close_PropEsts_groupSEs[i,seq(2,8,2)])
}
abline(v=4.5,lty=2,col="grey")

dev.off()	


## barplot
pdf("barplots/closeSingle_compEsts_viaSingleCell.pdf", w=50,h=6)
palette(pal)
par(mar = c(9,6,2,2), cex.axis=1.5,cex.lab=2)
ooClose= order(rse_gene_close$SampleLabel)
bp=barplot(t(close_PropEsts_Scaled[ooClose,]), col = pal,
	ylim = c(0,1.3),xaxt = "n", yaxt= "n",ylab="Class Proportion")

g = split(bp,rse_gene_close$SampleLabel[ooClose])
axis(1,at=sapply(g,mean), font=2, labels = gsub(" day ", "", names(g)),las=3)
abline(v = sapply(g,min)[-1]-0.5,lwd=4,col="black")
legend("top", colnames(close_PropEsts_Scaled), pch = 15, 
	col = 1:10,bg="white", nc = 4,cex=1.5,pt.cex=2)
axis(2, at=seq(0,1,by=0.25))
dev.off()

################
## camp fetal ##

# load data
load("/dcl01/ajaffe/data/lab/singleCell/camp_fetal/rse_gene_Camp_Fetal_Neurons_human_n734.Rdata")

## add pheno
pheno = read.delim("/dcl01/ajaffe/data/lab/singleCell/camp_fetal/camp_fetal_SraRunTable.txt",as.is=TRUE)
pheno = pheno[match(colnames(rse_gene), pheno$Run_s),]
colData(rse_gene) = cbind(colData(rse_gene), pheno)
rse_geneCamp = rse_gene
rse_geneCamp$Day = as.numeric(ss(rse_geneCamp$Stage_s, " "))
rse_geneCamp$Day[grep("weeks", rse_geneCamp$Stage_s)] = 7*rse_geneCamp$Day[grep("weeks", rse_geneCamp$Stage_s)] 
rse_geneCamp$Day[is.na(rse_geneCamp$Day)] = 12*7 # geo says 12 week pcw

yExprs_camp = log2(getRPKM(rse_geneCamp,"Length")+1)
yExprs_camp_Scaled = scale(yExprs_camp[rownames(coefEsts),])

camp_PropEsts= minfi:::projectCellType(yExprs_camp_Scaled, coefEsts)
camp_PropEsts_Scaled = prop.table(camp_PropEsts,1)


## camp et al data
gIndexes_camp = splitit(rse_geneCamp$Day)
camp_PropEsts_groupMeans = sapply(gIndexes_camp, 
	function(ii) colMeans(camp_PropEsts_Scaled[ii,]))
camp_PropEsts_groupSEs = sapply(gIndexes_camp, 
	function(ii) apply(camp_PropEsts_Scaled[ii,],2,function(x) sd(x)/sqrt(length(x))))

pdf("lineplots/camp_singleCell_cellTypeDecon.pdf",h=5,w=5)
palette(pal)
par(mar=c(10,6,2,2),cex.axis=1.2,cex.lab=1.3)
plot(camp_PropEsts_groupMeans[1,], type="b",xaxt="n",
	lty=1, pch=19,lwd=2, ylim = c(0,1),col=1,xlab="",
	ylab = "Class Proportion",cex=1.3)
axis(1,at=seq(along=gIndexes_camp),
	paste0(colnames(camp_PropEsts_groupMeans), "_PCD"), las=3,cex.axis=1.5)
segments(x0=seq(along=gIndexes_camp), 
	x1=seq(along=gIndexes_camp),col=1,lwd=2, 
	y0=camp_PropEsts_groupMeans[1,] - 2*camp_PropEsts_groupSEs[1,], 
	y1=camp_PropEsts_groupMeans[1,] + 2*camp_PropEsts_groupSEs[1,])

for(i in 2:10) {
	lines(camp_PropEsts_groupMeans[i,], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=seq(along=gIndexes_camp), 
		x1=seq(along=gIndexes_camp), col=i,lwd=2, 
	y0=camp_PropEsts_groupMeans[i,] - 2*camp_PropEsts_groupSEs[i,], 
	y1=camp_PropEsts_groupMeans[i,] + 2*camp_PropEsts_groupSEs[i,])
}
dev.off()	

pdf("barplots/campFetal_compEsts_viaSingleCell.pdf", w=50,h=6)
palette(pal)
par(mar = c(9,6,2,2), cex.axis=1.5,cex.lab=2)
ooCamp= order(rse_geneCamp$Day)
bp=barplot(t(camp_PropEsts_Scaled[ooCamp,]), col = pal,
	ylim = c(0,1.3),xaxt = "n", yaxt= "n",ylab="Class Proportion")

g = split(bp,rse_geneCamp$Day[ooCamp])
axis(1,at=sapply(g,mean), font=2, labels = gsub(" day ", "", names(g)),las=3)
abline(v = sapply(g,min)[-1]-0.5,lwd=4,col="black")
legend("top", colnames(camp_PropEsts_Scaled), pch = 15, 
	col = 1:10,bg="white", nc = 4,cex=1.5,pt.cex=2)
axis(2, at=seq(0,1,by=0.25))
dev.off()

#####################
## organoids ########

load("/dcl01/ajaffe/data/lab/singleCell/pasca_organoid/rse_gene_Pasca_Organoids_human_n730.Rdata")
rse_genePasca = rse_gene
yExprs_pasca = log2(getRPKM(rse_genePasca,"Length")+1)
yExprs_pasca_Scaled = scale(yExprs_pasca[rownames(coefEsts),])

pasca_PropEsts= minfi:::projectCellType(yExprs_pasca_Scaled,coefEsts)
pasca_PropEsts_Scaled = prop.table(pasca_PropEsts,1)

## find GEO ID
pdPasca_sra = read.delim("/dcl01/ajaffe/data/lab/singleCell/pasca_organoid/pasca_SraRunTable.txt",
	as.is=TRUE)
rse_genePasca$GSMID = pdPasca_sra$Sample_Name[match(colnames(rse_genePasca), pdPasca_sra$Run)]

## get GEO data
pascaData =GEOquery::getGEO("GSE99951")
pdPasca_geo = pData(pascaData[[1]])
pdPasca_geo = pdPasca_geo[rse_genePasca$GSMID,]

## seq type
rse_genePasca$SeqType = ifelse(grepl("bulk", pdPasca_geo$title), "Bulk", "Single")

## age
rse_genePasca$Age = as.numeric(ss(
	as.character(pdPasca_geo$characteristics_ch1.3), " ", 3))
rse_genePasca$Age[rse_genePasca$SeqType=="Bulk"] =  as.numeric(ss(ss(as.character(
	pdPasca_geo$title[rse_genePasca$SeqType=="Bulk"]), "_", 4), "b| "))
	
## cell type
rse_genePasca$CellType = ss(as.character(pdPasca_geo$characteristics_ch1.5), " ", 3)
rse_genePasca$CellType[pdPasca_geo$`immunopanning antibody:ch1` == "HepaCAM (astrocyte)"] = "glia"
rse_genePasca$CellType[pdPasca_geo$`immunopanning antibody:ch1` == "Thy1 (neuron)"] = "neuron"

# selection
rse_genePasca$Selection = ss(as.character(pdPasca_geo$characteristics_ch1.4), " ", 3)
rse_genePasca$Selection[pdPasca_geo$`immunopanning antibody:ch1` == "HepaCAM (astrocyte)"] = "HEPACAM"
rse_genePasca$Selection[pdPasca_geo$`immunopanning antibody:ch1` == "Thy1 (neuron)"] = "Thy1"

## split into single and bulk
rse_genePasca_single = rse_genePasca[,rse_genePasca$SeqType == "Single"]
rse_genePasca_bulk = rse_genePasca[,rse_genePasca$SeqType == "Bulk"]
pasca_PropEsts_Scaled_single = pasca_PropEsts_Scaled[rse_genePasca$SeqType == "Single",]
pasca_PropEsts_Scaled_bulk= pasca_PropEsts_Scaled[rse_genePasca$SeqType == "Bulk",]
pasca_PropEsts_single = pasca_PropEsts[rse_genePasca$SeqType == "Single",]
pasca_PropEsts_bulk = pasca_PropEsts[rse_genePasca$SeqType == "Bulk",]

## do single first
rse_genePasca_single$SampleLabel = paste0(substr(rse_genePasca_single$CellType, 1,3) , "_d", 
	rse_genePasca_single$Age, "_", gsub("HEPACAM", "HEP",gsub("Unselected", "Uns", rse_genePasca_single$Selection)))

## drop weird samples
keepIndex_single = which(! rse_genePasca_single$SampleLabel %in% c("gli_d450_Uns", "neu_d175_HEP"))
rse_genePasca_single = rse_genePasca_single[,keepIndex_single]
pasca_PropEsts_Scaled_single = pasca_PropEsts_Scaled_single[keepIndex_single,]
pasca_PropEsts_single = pasca_PropEsts_single[keepIndex_single,]

## relevel
rse_genePasca_single$SampleLabel = factor(rse_genePasca_single$SampleLabel,
	levels = names(table(rse_genePasca_single$SampleLabel))[c(1,3,2,5,4,6,7:11)])

## pasca data et al data
gIndexes_pasca_single = splitit(rse_genePasca_single$SampleLabel)
pasca_PropEsts_groupMeans_



#################
## load Tekin
library(GEOquery)

# load data
load("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/otherData/Tekin_nature/preprocessed_human/rse_gene_tekin_human_n157.Rdata")
rse_gene_tekin = rse_gene

## geo phenotypes
geoData = getGEO("GSE111831")
geo_pheno = do.call("rbind", lapply(geoData, function(x) pData(x)[,c("geo_accession","title")]))
rse_gene_tekin$geo_title = as.character(geo_pheno$title[match(rse_gene_tekin$Sample_Name, 
					geo_pheno$geo_accession)])
rse_gene_tekin$geo_label = ss(rse_gene_tekin$geo_title, "_rep")


type = c("human neuronal cells sorted from their 3D co-culture with human primary astrocytes",
		"human neuronal cells sorted from their 3D co-culture with differentiated human astrocytic cells",
		"human neuronal cells sorted from 3D only human neuronal cell culture",
		"human neuronal cells and mouse astrocytes",
		"human neuronal cells")
keepInd = which(rse_gene_tekin$strain_cell_type %in% type &
				rse_gene_tekin$culture_condition %in% c("2D","3D") &
				rse_gene_tekin$age %in% c("1wk","5wk"))

# keepInd = c(3:7,9:11,15:17,57:65,69:77,81:83)		
# keepInd = which(rse_gene_tekin$culture_condition %in% c("2D","3D") &
				# rse_gene_tekin$age %in% c("1wk","5wk") &
				# grepl("human neuronal cells", rse_gene_tekin$strain_cell_type))
				
rse_gene_tekin = rse_gene_tekin[,keepInd]

rse_gene_tekin$cells = ifelse(grepl("astro", rse_gene_tekin$strain_cell_type), "neur+ast", "neur")
rse_gene_tekin$SampleLabel =paste0(rse_gene_tekin$age, "_", rse_gene_tekin$culture_condition, "_", rse_gene_tekin$cells)

dropInd = which(rse_gene_tekin$SampleLabel == "5wk_3D_neur+ast" & 
			rse_gene_tekin$strain_cell_type == "human neuronal cells and mouse astrocytes")
rse_gene_tekin = rse_gene_tekin[,-dropInd]
table(pd$strain_cell_type, pd$SampleLabel)

pd = colData(rse_gene_tekin)
yExprs_tekin = log2(getRPKM(rse_gene_tekin,length_var="Length")+1)


# project
load("singleCell_iPSC_quake_coefEsts_calibration_Zscale.rda", verbose=TRUE)
yExprs_tekin_Z = scale(yExprs_tekin[rownames(coefEsts),])

propEsts_tekin = minfi:::projectCellType(yExprs_tekin_Z,coefEsts)
tekin_PropEsts_Scaled = prop.table(propEsts_tekin,1)

## line plot
gIndexes_tekin = splitit(rse_gene_tekin$SampleLabel)
tekin_PropEsts_groupMeans = sapply(gIndexes_tekin, 
	function(ii) colMeans(tekin_PropEsts_Scaled[ii,]))
tekin_PropEsts_groupSEs = sapply(gIndexes_tekin, 
	function(ii) apply(tekin_PropEsts_Scaled[ii,],2,function(x) sd(x)/sqrt(length(x))))


pdf("lineplots/tekin_singleCell_cellTypeDecon.pdf",h=5,w=5)
palette(pal)
par(mar=c(10,6,2,2),cex.axis=1.2,cex.lab=1.3)
plot(tekin_PropEsts_groupMeans[1,], type="p",xaxt="n",
	lty=1, pch=19,lwd=2, ylim = c(0,1),col=1,xlab="",
	ylab = "Class Proportion",cex=1.3)
## 1-2 3-4 5-6 7-8
segments(x0=seq(1,7,by=2), x1=seq(2,8,by=2), col=1,lwd=2, 
	y0=tekin_PropEsts_groupMeans[1,seq(1,7,by=2)], 
	y1=tekin_PropEsts_groupMeans[1,seq(2,8,by=2)])
	
axis(1,at=seq(along=gIndexes_tekin),
	colnames(tekin_PropEsts_groupMeans), las=3, cex.axis=1.2)
segments(x0=seq(along=gIndexes_tekin), 
	x1=seq(along=gIndexes_tekin),col=1,lwd=2, 
	y0=tekin_PropEsts_groupMeans[1,] - 2*tekin_PropEsts_groupSEs[1,], 
	y1=tekin_PropEsts_groupMeans[1,] + 2*tekin_PropEsts_groupSEs[1,])

for(i in 2:10) {
	lines(tekin_PropEsts_groupMeans[i,], type="p",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	## 1-2 3-4 5-6 7-8
	segments(x0=seq(1,7,by=2), x1=seq(2,8,by=2), col=i,lwd=2, 
		y0=tekin_PropEsts_groupMeans[i,seq(1,7,by=2)], 
		y1=tekin_PropEsts_groupMeans[i,seq(2,8,by=2)])
	segments(x0=seq(along=gIndexes_tekin), 
		x1=seq(along=gIndexes_tekin), col=i,lwd=2, 
		y0=tekin_PropEsts_groupMeans[i,] - 2*tekin_PropEsts_groupSEs[i,], 
		y1=tekin_PropEsts_groupMeans[i,] + 2*tekin_PropEsts_groupSEs[i,])
}
dev.off()	









boxplot(propEsts$Neuron ~ rse_gene$age_days)
boxplot(propEsts$Neuron ~ rse_gene$culture_condition)
plot(propEsts$Neuron ~ rse_gene$age_days)

par(mar=c(5,30,1,1))
boxplot(propEsts$iPSC ~ rse_gene$geo_label,horizontal=TRUE,las=1)
boxplot(propEsts$NPC ~ rse_gene$geo_label,horizontal=TRUE,las=1)
boxplot(propEsts$Fetal_quiescent ~ rse_gene$geo_label,horizontal=TRUE,las=1)
boxplot(propEsts$Neuron ~ rse_gene$geo_label,horizontal=TRUE,las=1)

boxplot(propEsts$Astrocytes ~ rse_gene$geo_label,horizontal=TRUE,las=1)
boxplot(propEsts$Oligodendrocytes ~ rse_gene$geo_label,horizontal=TRUE,las=1)
boxplot(propEsts$Endothelial ~ rse_gene$geo_label,horizontal=TRUE,las=1)

######################
## plots for paper  ##
######################

groups = c("2D_1wk_iN","2D_1wk_iN_ma", "2D_5wk_iN", "2D_5wk_iN_ma",
	"iN_5wk_3D_100ul_7.36mg/ml_Mtrgl_10mpml",
	"iN_from_5wk_3D_iN_astrocytic-cell_100ul_7.36mg/ml_Mtrgl_20mpml",
	"iN_from_5wk_3D_iN_human-primary-astro_100ul_7.36mg/ml_Mtrgl_20mpml",
	



