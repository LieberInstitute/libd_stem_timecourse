###
library(SummarizedExperiment)
library(minfi)
library(recount)

pd = read.csv("metadata_56cells.csv",as.is=TRUE)

load("rse_gene_Bardy_singleCell_n56.Rdata")
rse_gene = rse_gene[,pd$cell.number]

## do decon
load("../../deconvolution/cell_type/singleCell_iPSC_quake_coefEsts_calibration_Zscale.rda")
geneExprs= log2(getRPKM(rse_gene, "Length")+1)
geneExprs_scale = scale(geneExprs[rownames(coefEsts),])

cellPropEsts = minfi:::projectCellType(geneExprs_scale,coefEsts)
cellPropEsts = as.data.frame(cellPropEsts)
pd = cbind(pd, cellPropEsts)

write.csv(pd, file = "metadata_56cells_withRnaComp.csv", row.names=FALSE)

pdf("bardy_ePhys_comp_checks.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2, cex.main=2)
boxplot(pd$Neuron ~ pd$FTC_simple, outline=FALSE,
	ylab = "Neuronal RNA Fraction", xlab="AP Type")
points(pd$Neuron ~ jitter(as.numeric(factor(pd$FTC_simple))),
	pch = 21, bg = "grey", cex=2)

boxplot(pd$Fetal_quiescent ~ pd$FTC_simple,outline=FALSE,
	ylab = "Fetal Quiescent RNA Fraction", 
	xlab="AP Type")
points(pd$Fetal_quiescent ~ jitter(as.numeric(factor(pd$FTC_simple))),
	pch = 21, bg = "grey", cex=2)

boxplot(pd$NPC ~ pd$FTC_simple,outline=FALSE,
	ylab = "NPC RNA Fraction", xlab="AP Type")
points(pd$NPC ~ jitter(as.numeric(factor(pd$FTC_simple))),
	pch = 21, bg = "grey", cex=2)
boxplot(pd$Endothelial ~ pd$FTC_simple,outline=FALSE,
	ylab = "Endothelial RNA Fraction", xlab="AP Type")
points(pd$Endothelial ~ jitter(as.numeric(factor(pd$FTC_simple))),
	pch = 21, bg = "grey", cex=2)
dev.off()

boxplot(pd$NPC ~ pd$FTC_simple)
plot(pd$Neuron ~ pd$FTC_simple)

summary(lm(pd$Neuron ~ pd$FTC_simple, subset=pd$FTC_simple!=0))
summary(lm(pd$Neuron ~ pd$FTC, subset=pd$FTC_simple!=0))
summary(lm(pd$Fetal_quiescent ~ pd$FTC_simple, subset=pd$FTC_simple!=0))
summary(lm(pd$Fetal_replicating ~ pd$FTC_simple, subset=pd$FTC_simple!=0))
summary(lm(pd$NPC ~ pd$FTC_simple))
