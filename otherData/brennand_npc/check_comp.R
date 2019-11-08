##
library(jaffelab)
library(minfi)
library(SummarizedExperiment)
library(recount)
library(limma)
library(edgeR)
library(RColorBrewer)

## read in gene
geneCounts = read.csv("geneCounts.csv", as.is=TRUE)
# pd = read.csv("info.csv",as.is=TRUE)
pd = read.csv("COS RNAseq Master Sheet - Co-variates for analysis.csv",as.is=TRUE)
geneMap = read.delim("ENSEMBLv70_gene_info.tsv",as.is=TRUE,row.names=1)

## fix name
pd$SampleID = paste0("X", gsub("-", ".", pd$Sample.Name))
colnames(geneCounts) = gsub("N.2$", "NR", colnames(geneCounts))
colnames(geneCounts) = gsub("F.2$", "FR", colnames(geneCounts))

pd$Cell.Type = ifelse(pd$Cell.Type == "NPC", "NPC", "Neuron")
pd$Cell.Type = factor(pd$Cell.Type , levels =c("NPC", "Neuron"))

mm = match(pd$SampleID,colnames(geneCounts))
table(table(mm)) # fixed

# reorder
geneCounts = geneCounts[,pd$SampleID]
identical(rownames(geneCounts), rownames(geneMap)) # TRUE
identical(colnames(geneCounts), pd$SampleID) # TRUE

#############
## drop #####
#############

keepIndex = which(pd$Exclusion.Criteria == "")
pd = pd[keepIndex,]
geneCounts = geneCounts[,keepIndex]

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
 
pdf("../../deconvolution/cell_type/lineplots/brennand_stemcell_timecourse_cellTypeDecon.pdf",h=5,w=5)
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

## merge
pd = cbind(pd, cellPropEsts)

## fix donor ID
pd$Donor = ss(pd$Sample.Name, "-")

#####################
## do DE modeling ###
#####################

## filter low expression
gIndex = rowMeans(geneRpkm) > 0.2
geneCounts = geneCounts[gIndex,]
geneMap = geneMap[gIndex,]

## by cell type
cIndexes = splitit(pd$Cell.Type)

dgeList= lapply(cIndexes, function(ii) {
	dge = DGEList(counts = geneCounts[,ii], genes = geneMap)
	calcNormFactors(dge)
})

modList = lapply(cIndexes, function(ii) {
	mod = model.matrix(~ Dx + Sex + NPC + Fetal_quiescent + Neurons, data=pd[ii,])
})

voomList = mapply(function(d,m) voom(d, m, plot = FALSE), dgeList,modList,SIMPLIFY=FALSE)

# corObjList = mclapply(1:2,function(i) {
	# duplicateCorrelation(voomList[[i]]$E, modList[[i]], block=pd$Donor[cIndexes[[i]]])
# },mc.cores=2)
# save(corObjList, file = "duplCorr_geneLevel_byCellType.rda")
load("duplCorr_geneLevel_byCellType.rda")

fitList = mclapply(1:2,function(i) {
	geneFit = lmFit(voomList[[i]]$E, modList[[i]],
		block=pd$Donor[cIndexes[[i]]],
		correlation = corObjList[[i]]$consensus.correlation)
	eBayes(geneFit)
},mc.cores=2)
names(fitList) = names(cIndexes)

### get stats
szStatList = lapply(fitList, topTable, coef = 2, n = nrow(geneCounts), sort="none")
szStats = do.call("cbind", szStatList)
szStats$Symbol = geneMap$geneName

plot(szStats$NPC.t, szStats$Neuron.t)
colSums(szStats[,grep("adj", colnames(szStats))] < 0.1)

######################
## interaction model #
######################

modListInt = modList
modListInt$NPC = model.matrix(~ Dx * NPC + Sex +  Fetal_quiescent + Neurons, data=pd[cIndexes$NPC,])
modListInt$Neuron = model.matrix(~ Dx * Neurons  + Sex +  NPC + Fetal_quiescent, data=pd[cIndexes$Neuron,])

fitListInt = mclapply(1:2,function(i) {
	geneFit = lmFit(voomList[[i]]$E, modListInt[[i]],
		block=pd$Donor[cIndexes[[i]]],
		correlation = corObjList[[i]]$consensus.correlation)
	eBayes(geneFit)
},mc.cores=2)
names(fitListInt) = names(cIndexes)

## extract
szStatMainList = lapply(fitListInt, topTable, coef = 2, 
	n = nrow(geneCounts), sort="none")
szStatsMain = do.call("cbind", szStatMainList)
colnames(szStatsMain) = paste0(colnames(szStatsMain), "_main")
szStatIntList = lapply(fitListInt, topTable, coef = 7, 
	n = nrow(geneCounts), sort="none")
szStatsInt = do.call("cbind", szStatIntList)
colnames(szStatsInt) = paste0(colnames(szStatsInt), "_int")

colSums(szStatsMain[,grep("adj", colnames(szStatsMain))] < 0.05)
colSums(szStatsInt[,grep("adj", colnames(szStatsMain))] < 0.05)

## merge
szStatsEx = cbind(szStatsMain, szStatsInt)
szStatsEx$Symbol = geneMap$geneName
plot(szStatsEx$NPC.t_main, szStats$NPC.t)

sigStats = szStatsEx[rowSums(szStatsEx[,grep("adj", colnames(szStatsEx))] < 0.05) > 0,]
sigStats = sigStats[order(sigStats$NPC.P.Value_main),!grepl("Ave|B", colnames(sigStats))]
write.csv(sigStats, file = "hoffman_cellModeling_sigStats.csv")

origStats = szStats[rownames(sigStats),]

#######################
## show some plots?? ##
pd$DonorFac = factor(pd$Donor, levels = unique(pd$Donor[order(pd$Dx,pd$Donor)]))
geneExprs = as.matrix(geneExprs)
npcInd = pd$Cell.Type == "NPC"

pd$CellQuad = NA
pd$CellQuad[pd$Cell.Type == "NPC"] = cut(pd$NPC[pd$Cell.Type == "NPC"], 
	breaks= quantile(pd$NPC[pd$Cell.Type == "NPC"]), include=TRUE,labels=FALSE)
pd$CellQuad[pd$Cell.Type == "Neuron"] = cut(pd$Neurons[pd$Cell.Type == "Neuron"], 
	breaks= quantile(pd$Neurons[pd$Cell.Type == "Neuron"]), include=TRUE,labels=FALSE)

pd$DxQuad = paste0(pd$CellQuad, ":", pd$Dx)
pdf("example_boxplots_of_interaction.pdf",width=10)
par(mar=c(5,6,4,2), cex.axis=1.8,cex.lab=1.8,cex.main=1.8)
palette(c("blue", "red"))
for(i in 1:nrow(sigStats)) {
	boxplot(geneExprs[rownames(sigStats)[i],npcInd] ~ pd$DxQuad[npcInd],
		main = paste0("NPC - ", sigStats$Symbol[i]), ylab= "log2(RPKM+1)",
		xlab = "NPC RNA Fraction Quantile & Dx", # col=rep(1:2,times=4),
		outline=FALSE, ylim = range(geneExprs[rownames(sigStats)[i],npcInd]))
	points(geneExprs[rownames(sigStats)[i],npcInd] ~ 
		jitter(as.numeric(factor(pd$DxQuad)[npcInd]), amount=0.15), 
		pch = 21, bg=factor(ss(pd$DxQuad[npcInd], ":", 2)))
	abline(v = c(2.5,4.5,6.5),lty=2,cex=1.8)
	legend("topleft", paste0(c("Main", "Int"),
		" p=", signif(sigStats[i,c(3,11)],3)),nc=1,cex=1.8)
}
dev.off()

plot(geneExprs[rownames(sigStats)[1],npcInd] ~ pd$NPC[npcInd])
cor.test(geneExprs[rownames(sigStats)[1],npcInd] ,  pd$NPC[npcInd])


		
## model matrix 
geneFit= lmFit(vGene, mod)
geneFit= eBayes(geneFit)
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