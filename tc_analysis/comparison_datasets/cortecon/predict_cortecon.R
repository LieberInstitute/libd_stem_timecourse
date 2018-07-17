###
library(minfi)
library(jaffelab)
library(SummarizedExperiment)
library(genefilter)
library(edgeR)
library(preprocessCore)

## load cortecon data ## hg38
load("/dcs01/ajaffe/LIBD_StemCell/Data/processed_rawCounts_Cortecon.rda") ## for Day, other pd
pdCort2 = pd[,1:8]
load("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/cortecon/rawCounts_cortecon_hg38_n24.rda")
pdCort = metrics
stopifnot(identical(rownames(pdCort),pdCort2$SRR))
pdCort = cbind(pdCort,pdCort2)
geneCountsCort = geneCounts

## load  ## hg38
MAINDIR = '/dcl01/lieber/ajaffe/lab/libd_stem_timecourse'
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n146.rda"))
geneCounts = assays(rse_gene)$counts
pd = colData(rse_gene)

## don't use technical replicates
bioInd = which(pd$Class == "Naked genomes")
pd = pd[bioInd,]
geneCounts = geneCounts[,bioInd]   # n=128
## don't use RENEW controls or NEURONS_ALONE in analyses
dropInd = which(pd$CONDITION %in% c("RENEW","NEURONS_ALONE"))
pd = pd[-dropInd,]
geneCounts = geneCounts[,-dropInd]   # n=106


## specify study
study = rep(c("Cort","LIBD"), 
	times = c(nrow(pdCort), nrow(pd)))
sIndexes = splitit(study)

## normalize together
counts = cbind(geneCountsCort,geneCounts)
countsNorm = normalize.quantiles(as.matrix(counts))
rownames(countsNorm) = rownames(geneCounts)

### split up
yCort = as.matrix(log2(countsNorm[,study=="Cort"]+1))

cIndexes = splitit(pdCort$Day)
tstatList <- lapply(cIndexes, function(i) {
	x <- rep(0,ncol(yCort))
	x[i] <- 1
	return(rowttests(yCort, factor(x)))
})

# pick cell type specific genes
numGenes = 5
geneList <- lapply(tstatList, function(x) {
	y <- x[x[,"p.value"] < 1e-6,]
	yUp <- y[order(y[,"dm"], decreasing=TRUE),]
	yDown <- y[order(y[,"dm"], decreasing=FALSE),]
	c(rownames(yUp)[1:numGenes], rownames(yDown)[1:numGenes])
})

compGenes = unique(unlist(geneList))
yComp = yCort[compGenes,]
mapComp = geneMap[compGenes,]

## retain means for batch effects check
pdCort$sampleLog2Means = colMeans(yComp)

## modeling
pdCort$DayFactor= factor(paste0("D",pdCort$Day),
	levels = paste0("D", sort(unique(pdCort$Day))))

form <- as.formula(sprintf("y ~ %s - 1", 
	paste(levels(pdCort$DayFactor), 	collapse="+")))
phenoDF <- as.data.frame(model.matrix(~pdCort$DayFactor-1))
colnames(phenoDF) <- gsub(" ", "_", sub("^pdCort\\$DayFactor", "", colnames(phenoDF)))

### estimate
coefEsts <- minfi:::validationCellType(Y = yComp, 
	pheno = phenoDF, modelFix = form)$coefEsts
save(coefEsts, file = "cortecon_stemCell_dayEstimates_QN.rda")

###################
# load stem

yRpkm = log2(countsNorm[,study=="LIBD"]+1)
dayEstimates = minfi:::projectCellType(
	yRpkm[rownames(coefEsts),], coefEsts)
dayEstimates = as.data.frame(dayEstimates)
save(dayEstimates, file="LIBD_dayEstimates_fromCortecon.rda")

pd$sampleMeans = colMeans(yRpkm)
plot(c(pdCort$sampleLog2Means, pd$sampleMeans))

pd$day = as.numeric(pd$DAY)
pd$Donor[pd$Donor == "NA"] = NA
pd$Donor = factor(pd$Donor)

library(RColorBrewer)
pal = c(colorRampPalette(brewer.pal(10,"Paired")[1:2])(4), 
		colorRampPalette(brewer.pal(10,"Paired")[3:4])(3), 
		brewer.pal(10,"Paired")[9:10], brewer.pal(10,"Paired")[5:8])
pd$lineCol = pal[as.numeric(factor(pd$LINE))]
donorInd = c(rep(4,4),rep(7,3),9,9,11,11,13,13)
pd$donorCol = pal[donorInd[as.numeric(factor(pd$LINE))]]


pdf("cortecon_estimates_QN.pdf",w=8,h=6)
par(mar=c(10,6,2,2), cex.axis=2, cex.lab=2, cex.main=2)
for(i in 1:ncol(dayEstimates)) {
	boxplot(dayEstimates[,i] ~ pd$day, ylim = c(0,1),
		xlab="LIBD Stem Day",outline=FALSE,
		ylab = "Proportion")
	points(dayEstimates[,i] ~ jitter(as.numeric(factor(pd$day)),
		amount=0.15), pch = 21, bg=pd$donorCol, cex=2)
	legend("top",paste("Cortecon", colnames(dayEstimates)[i]), cex=2, text.font=2, bty = "n")
	if (i==1) { legend("topright", paste0("Donor ",levels(pd$Donor)),
		col = pal[c(4,7,9,11,13)], pch = 15, cex=1.1) }
}
dev.off()

pd$COND = pd$CONDITION
pd$COND[grep("ACC_DORS",pd$COND)] = "ACC_DORSAL"
pd$COND = as.factor(pd$COND)
pd$COND = factor(pd$COND,levels(pd$COND)[c(1,3,4,2)])  ## put levels in order

pdf("cortecon_estimates_condition.pdf",w=8, h=6)
par(mar=c(10,6,2,2), cex.axis=1.5, cex.lab=2, cex.main=2)
for (i in 1:ncol(dayEstimates)) {
	boxplot(dayEstimates[,i] ~ pd$COND,outline=FALSE, ylim=c(0,1), xaxt="n")
	points(dayEstimates[,i] ~ jitter(as.numeric(pd$COND),amount=0.15), pch = 21, bg=pd$donorCol, cex=2)
	axis(1, at=1:4, labels = FALSE)
	text(1:4, par("usr")[3]-.07,
     srt = 40, adj= 1, xpd = TRUE,
     labels = levels(pd$COND), cex=1.25)
	legend("top",paste("Cortecon", colnames(dayEstimates)[i]), cex=2, text.font=2, bty = "n")
	if (i==1) { legend("topright", paste0("Donor ",levels(pd$Donor)),
		col = pal[c(4,7,9,11,13)], pch = 15, cex=1.1) }
}
dev.off()




pdf("cortecon_estimates_day.pdf",w=9)
palette(brewer.pal(8,"Dark2"))
par(mar=c(5,6,4,2), cex.axis=1, cex.lab=1.8, cex.main=1.8)
for (i in 1:ncol(dayEstimates)) {
boxplot(dayEstimates[,i] ~ pd$COND,outline=FALSE, ylim=c(0,1), main=colnames(dayEstimates)[i])
	points(dayEstimates[,i] ~ jitter(as.numeric(pd$COND),amount=0.15), pch = 21, bg=factor(pd$Donor))
}
dev.off()

plot(dayEstimates$D7 ~ pd$DAY)
plot(dayEstimates$D63 ~ pd$DAY)





