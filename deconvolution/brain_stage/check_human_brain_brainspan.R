#####

library(SummarizedExperiment)
library(recount)
library(jaffelab)
library(RColorBrewer)
library(minfi)

MAINDIR = "/dcl01/lieber/ajaffe/lab/libd_stem_timecourse"

# ## human data from span
# load("/dcl01/lieber/ajaffe/lab/brainspan_analysis/rawCounts_full_brainspan_n607.rda", verbose=TRUE)
# load("/dcl01/lieber/ajaffe/lab/brainspan_analysis/full.brainspan_final.phenotype.rda",verbose=TRUE)
# pdSpan = cbind(fin_pdSpan, metrics[fin_pdSpan$lab,])
# pdSpan$wig = NULL
# geneCounts = geneCounts[,pdSpan$lab]
# rownames(pdSpan) = pdSpan$lab
# rm(exonCounts, exonMap, jCounts, jMap)

# ## Create gene,exon RangedSummarizedExperiment objects
# gr_genes <- GRanges(seqnames = geneMap$Chr,
    # IRanges(geneMap$Start, geneMap$End), strand = geneMap$Strand)
# names(gr_genes) <- rownames(geneMap)
# mcols(gr_genes) <- DataFrame(geneMap[, - which(colnames(geneMap) %in%
    # c('Chr', 'Start', 'End', 'Strand'))])
# rse_geneSpan <- SummarizedExperiment(assays = list('counts' = geneCounts),
    # rowRanges = gr_genes, colData = pdSpan)
	
# ### using AgeCat
# rse_geneSpan$AgeCat = "Adult"
# rse_geneSpan$AgeCat[rse_geneSpan$Age<20] = "Teens"
# rse_geneSpan$AgeCat[rse_geneSpan$Age<10] = "Child"
# rse_geneSpan$AgeCat[rse_geneSpan$Age<1] = "Infant"
# rse_geneSpan$AgeCat[rse_geneSpan$Age<0] = "Late Fetal"
# rse_geneSpan$AgeCat[rse_geneSpan$Age<(-.25)] = "Mid Fetal"
# rse_geneSpan$AgeCat[rse_geneSpan$Age<(-.5)] = "Early Fetal"
# rse_geneSpan$AgeCat = as.factor(rse_geneSpan$AgeCat)
# rse_geneSpan$AgeCat = factor(rse_geneSpan$AgeCat,levels(rse_geneSpan$AgeCat)[c(3,6,5,4,2,7,1)])

# ## brainspan NCX
# rse_geneSpan$RegionGroup = rse_geneSpan$Regioncode 
# rse_geneSpan$RegionGroup[rse_geneSpan$Regioncode %in% 
				# c("A1C", "DFC", "IPC","ITC", "M1C", "MFC",
				# "OFC", "S1C", "STC", "V1C", "VFC")] = "NCX"

# save(rse_geneSpan, file = "brainspan_gene_rse.rda")
load("brainspan_gene_rse.rda")

geneRpkmSpan = getRPKM(rse_geneSpan,length_var="Length")
yRpkmSpan= log2(geneRpkmSpan+1)

################
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


############
## filter and merge

## just iPSC
rse_geneSC = rse_geneSC[,rse_geneSC$CONDITION == "iPSC"] #
yExprsSC = log2(getRPKM(rse_geneSC, "Length")+1)

# just neocortex
keepIndex = which(rse_geneSpan$RegionGroup == "NCX")
rse_geneSpan_filter = rse_geneSpan[,keepIndex]
yExprsSpan_filter = yRpkmSpan[,keepIndex]
		
rse_geneSpan_filter$AgeRegion = paste0(rse_geneSpan_filter$RegionGroup,"_", 
					gsub(" ", "", rse_geneSpan_filter$AgeCat))
rse_geneSpan_filter$AgeRegion = factor(rse_geneSpan_filter$AgeRegion,
	paste0("NCX_", 	gsub(" ", "", levels(rse_geneSpan$AgeCat))))		
	
## merge w/ scorecard
yExprs_Merge = cbind(yExprsSC, yExprsSpan_filter)	
group = c(as.character(rse_geneSC$CONDITION), 
	as.character(rse_geneSpan_filter$AgeRegion))
group = factor(group,levels =unique(group))

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
mergeMarkerExprs <- yExprs_Merge[trainingProbes, ]
mergeMarkerMeanExprs <- colMeans(mergeMarkerExprs)
names(mergeMarkerMeanExprs) <- names(tIndexes)

form <- as.formula(sprintf("y ~ %s - 1", paste(levels(group),collapse = "+")))
phenoDF <- as.data.frame(model.matrix(~group - 1))
colnames(phenoDF) <- sub("^group", "", colnames(phenoDF))

## try z-score
mergeMarkerExprsZ = scale(mergeMarkerExprs)
mergeMarkerMeanExprsZ = colMeans(mergeMarkerExprsZ)

## do calibration
coefEsts <- minfi:::validationCellType(Y = mergeMarkerExprsZ, 
	pheno = phenoDF, modelFix = form)$coefEsts
save(coefEsts, file = "iPSC_BrainSpan_coefEsts_calibration_Zscore.rda")
write.csv(coefEsts, file="iPSC_BrainSpan_coefEsts_calibration_Zscore.csv")

##################################
### project into brainseq ########
##################################
load(file = "iPSC_BrainSpan_coefEsts_calibration_Zscore.rda")

pal = brewer.pal(8,"Dark2")

### validate using BrainSeq
load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseGene_merged_n732.rda")
rse_geneDLPFC = rse_gene[,rse_gene$Dx == "Control" & 
	rse_gene$Race %in% c("AA","CAUC")]
rse_geneDLPFC$ageGroup = cut(rse_geneDLPFC$Age, c(-1,-0.45, 0,1,10,20,100))
levels(rse_geneDLPFC$ageGroup) = c("EarlyFetal","MidFetal",
	"Infant","Child","Teens","Adult")

rse_geneDLPFC = rse_geneDLPFC[,order(rse_geneDLPFC$Age)]
yExprsDLPFC = log2(getRPKM(rse_geneDLPFC,length_var="Length")+1)

yExprsDLPFC_Z = scale(yExprsDLPFC[rownames(coefEsts),])
dlpfcPropEsts = minfi:::projectCellType(yExprsDLPFC_Z,coefEsts)

# project
dlpfcPropEsts = minfi:::projectCellType(yExprsDLPFC_Z,coefEsts)

# summarize
gIndexes_dlpfc = splitit(rse_geneDLPFC$ageGroup)
dlpfcPropEsts_groupMeans = sapply(gIndexes_dlpfc, 
	function(ii) colMeans(dlpfcPropEsts[ii,]))
dlpfcPropEsts_groupSEs = sapply(gIndexes_dlpfc, 
	function(ii) apply(dlpfcPropEsts[ii,],2,function(x) sd(x)/sqrt(length(x))))

## make legend pdf
pdf("lineplots/legend_brainStageDecon.pdf")
plot(0,0, axes=FALSE,type="n")
legend("top", gsub("NCX_", "", rownames(dlpfcPropEsts_groupMeans)),
	ncol = 1, pch = 15, col = 1:8,pt.cex=1.5,cex=1.3)
dev.off()

## make line plot
pdf("lineplots/brainseq_phase1_dlpfc_brainStageDecon.pdf",h=5,w=5)
palette(pal)
par(mar=c(10,6,2,2),cex.axis=1.2,cex.lab=1.5)
plot(dlpfcPropEsts_groupMeans[1,], type="b",xaxt="n",
	lty=1, pch=19,lwd=2, ylim = c(0,1),col=1,xlab="",
	ylab = "Class Proportion",cex=1.3)
axis(1,at=1:6,colnames(dlpfcPropEsts_groupMeans), las=3,cex.axis=1.5)
segments(x0=1:7, x1=1:7, col=1,lwd=2, 
	y0=dlpfcPropEsts_groupMeans[1,] - 2*dlpfcPropEsts_groupSEs[1,], 
	y1=dlpfcPropEsts_groupMeans[1,] + 2*dlpfcPropEsts_groupSEs[1,])

for(i in 2:8) {
	lines(dlpfcPropEsts_groupMeans[i,], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=1:7, x1=1:7, col=i,lwd=2, 
	y0=dlpfcPropEsts_groupMeans[i,] - 2*dlpfcPropEsts_groupSEs[i,], 
	y1=dlpfcPropEsts_groupMeans[i,] + 2*dlpfcPropEsts_groupSEs[i,])
}
dev.off()


###################
## based on stem ##
###################
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n146.rda"))
## don't use technical replicates
bioInd = which(colData(rse_gene)$Class == "Naked genomes")
rse_gene = rse_gene[,bioInd]  # n=128
colData(rse_gene) = colData(rse_gene)
colData(rse_gene)$DAY = paste0("day_",colData(rse_gene)$DAY)
colData(rse_gene)$experiment = "lieber"
rse_geneTC = rse_gene

## clean CONDITION
rse_geneTC$COND = rse_geneTC$CONDITION
rse_geneTC$COND[grep("ACC_DORSAL",rse_geneTC$COND)] = "ACC_DORSAL"
rse_geneTC$COND[grep("PLUS",rse_geneTC$COND)] = "NEURONS+ASTROS"
rse_geneTC$COND = as.factor(rse_geneTC$COND)
rse_geneTC$COND = factor(rse_geneTC$COND,levels(rse_geneTC$COND)[c(5,1,4,6,2,3)])

## rpkm function
getRPKM = recount::getRPKM

geneRpkmTC = getRPKM(rse_geneTC, "Length")
yRpkmTC = log2(geneRpkmTC+1)
yRpkmTC_Z = scale(yRpkmTC[rownames(coefEsts),])

### fit model
stemPropEsts = minfi:::projectCellType(yRpkmTC_Z,coefEsts)
stemPropEsts_Prop = prop.table(stemPropEsts,1)

### plots
gIndexes_tc = splitit(rse_geneTC$COND)
stemPropEsts_groupMeans = sapply(gIndexes_tc, 
	function(ii) colMeans(stemPropEsts[ii,]))
stemPropEsts_groupSEs = sapply(gIndexes_tc, 
	function(ii) apply(stemPropEsts[ii,],2,function(x) sd(x)/sqrt(length(x))))

# for the paper
round(100*stemPropEsts_groupMeans,2)


pdf("lineplots/libd_stemcell_timecourse_brainStageDecon.pdf",h=5,w=5)
palette(pal)
par(mar=c(10,6,2,2),cex.axis=1.5,cex.lab=1.5)
plot(stemPropEsts_groupMeans[1,], type="b",xaxt="n",
	lty=1, pch=19,lwd=2, ylim = c(0,1),col=1,xlab="",
	ylab = "Class Proportion",cex=1.3)
axis(1,at=seq(along=gIndexes_tc),c("RENEW","ACCEL\nDORSAL", "NPC", "ROSETTE", 
	"NEURONS","NEURONS\n+ASTROS"), las=3,cex.axis=1.5)
segments(x0=1:7, x1=1:7, col=1,lwd=2, 
	y0=stemPropEsts_groupMeans[1,] - 2*stemPropEsts_groupSEs[1,], 
	y1=stemPropEsts_groupMeans[1,] + 2*stemPropEsts_groupSEs[1,])

for(i in 2:8) {
	lines(stemPropEsts_groupMeans[i,], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=1:7, x1=1:7, col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,] - 2*stemPropEsts_groupSEs[i,], 
	y1=stemPropEsts_groupMeans[i,] + 2*stemPropEsts_groupSEs[i,])
}
dev.off()

## barplot just of neurons
neuronIndex =grep("^NEURON", rse_geneTC$CONDITION )
pdf("barplots/LIBDstem_compEsts_justNeurons.pdf",w=4,h=5)
palette(pal)
par(mar=c(10,6,2,2),cex.axis=1.3,cex.lab=1.5)
bp = barplot(t(stemPropEsts_Prop[neuronIndex,]), col = pal,
	ylim = c(0,1),xaxt = "n", yaxt= "n",ylab="Class Proportion")
g = split(bp, rse_geneTC$COND[neuronIndex])
axis(1,at=sapply(g,mean)[5:6], font=1, labels = c("",""))
abline(v = sapply(g,min)[6]-0.5,lwd=2,col="black")
axis(2, at=seq(0,1,by=0.25))
dev.off()


## test for differences
modAdd = model.matrix(~ as.numeric(COND) + LINE + totalAssignedGene, 
	data=colData(rse_geneTC))
# > signif(t(apply(stemPropEsts, 2, function(y) summary(lm(y ~ modAdd - 1))$coef[2,])),3)             Estimate Std. Error t value Pr(>|t|)
# iPSC           -1.31e-01   6.71e-03 -19.600 3.77e-38
# NCX_EarlyFetal  4.40e-02   8.23e-03   5.350 4.70e-07
# NCX_MidFetal    5.50e-02   4.53e-03  12.200 2.97e-22
# NCX_LateFetal   1.65e-02   3.37e-03   4.910 3.13e-06
# NCX_Infant      7.95e-03   1.50e-03   5.300 5.71e-07
# NCX_Child      -6.98e-19   5.48e-18  -0.127 8.99e-01
# NCX_Teens      -6.45e-18   3.88e-18  -1.660 9.94e-02
# NCX_Adult       2.99e-02   2.29e-03  13.000 2.74e-24


mod = model.matrix(~ COND + LINE + totalAssignedGene, 
	data=colData(rse_geneTC))
statList = lapply(as.data.frame(stemPropEsts), function(y) summary(lm(y ~ mod - 1))$coef[2:6,])
tMat = t(sapply(statList, function(x) x[,3]))
pMat = t(sapply(statList, function(x) x[,4]))
colnames(tMat) = colnames(pMat) = gsub("modCOND", "", colnames(tMat))

## switches
switchStatList = lapply(1:(length(cIndexes)-1), function(i) {
	ii = unlist(cIndexes[i:(i+1)])
	cond = rse_geneTC$COND[ii]
	cond = factor(as.character(cond), 
		levels = levels(rse_geneTC$COND)[levels(rse_geneTC$COND) %in% cond])
	m = model.matrix(~ cond + LINE + totalAssignedGene, 
		data=colData(rse_geneTC[,ii]))
	t(sapply(as.data.frame(stemPropEsts[ii,]), 
		function(y) summary(lm(y ~ m - 1))$coef[2,]))
})	
pMatSwitch =  sapply(switchStatList, function(x) x[,4])
tMatSwitch =  sapply(switchStatList, function(x) x[,3])

colnames(pMatSwitch) = colnames(tMatSwitch) = sapply(1:(length(cIndexes)-1), function(i) {
	ii = unlist(cIndexes[i:(i+1)])
	cond = rse_geneTC$COND[ii]
	paste(levels(rse_geneTC$COND)[levels(rse_geneTC$COND) %in% cond], collapse=">")
})

round(tMatSwitch,2)
signif(pMatSwitch,3)





## stem
rse_geneTC$DAYNUM = as.numeric(ss(rse_geneTC$DAY, "_",2))
pdf("LIBDstem_compEsts.pdf",w=14,h=8)
palette(pal)
par(mar = c(14,6,2,2), cex.axis=1.5,cex.lab=2)
ooStem = order(rse_geneTC$COND, rse_geneTC$DAYNUM)
bp = barplot(t(stemPropEsts[ooStem,]), col = pal,
	ylim = c(0,1.3),xaxt = "n", yaxt= "n",ylab="Class Proportion")
g = split(bp, rse_geneTC$COND[ooStem])

axis(1,at=sapply(g,mean), font=2, labels = names(g),las=3)
abline(v = sapply(g,min)[-1]-0.5,lwd=2,col="black")
legend("top", colnames(stemPropEsts), pch = 15, 
	col = 1:8,bg="white", nc = 4,cex=1.5,pt.cex=2)
axis(2, at=seq(0,1,by=0.25))
dev.off()

## dlpfc
pdf("brainSeq_compEsts.pdf", w=14,h=6)
palette(pal)
par(mar = c(6,6,2,2), cex.axis=1.5,cex.lab=2)
ooDlpfc = order(rse_geneDLPFC$Age)
bp=barplot(t(dlpfcPropEsts[ooDlpfc,]), col = pal,
	ylim = c(0,1.3),xaxt = "n", yaxt= "n",ylab="Class Proportion")

g = split(bp, rse_geneDLPFC$ageGroup[ooDlpfc])
axis(1,at=sapply(g,mean), font=2, labels = names(g),las=3)
abline(v = sapply(g,min)[-1]-0.5,lwd=4,col="black")
legend("top", colnames(dlpfcPropEsts), pch = 15, 
	col = 1:8,bg="white", nc = 4,cex=1.5,pt.cex=2)
axis(2, at=seq(0,1,by=0.25))
dev.off()

################################
##### check other datasets #####

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
spanSingle_PropEsts= minfi:::projectCellType(yExprs_SpanSingle[trainingProbes,],coefEsts)


## brainspan single cell
pdf("brainSpanSingle_compEsts.pdf", w=30,h=6)
palette(pal)
par(mar = c(6,6,2,2), cex.axis=1.5,cex.lab=2)
ooSpanSingle = order(rse_geneSpanSingle$PCW)
bp=barplot(t(spanSingle_PropEsts[ooSpanSingle,]), col = pal,
	ylim = c(0,1.3),xaxt = "n", yaxt= "n",ylab="Class Proportion")

g = split(bp,rse_geneSpanSingle$PCW[ooSpanSingle])
axis(1,at=sapply(g,mean), font=2, labels = names(g),las=3)
abline(v = sapply(g,min)[-1]-0.5,lwd=4,col="black")
legend("top", colnames(spanSingle_PropEsts), pch = 15, 
	col = 1:8,bg="white", nc = 4,cex=1.5,pt.cex=2)
axis(2, at=seq(0,1,by=0.25))
dev.off()

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
close_PropEsts= minfi:::projectCellType(yExprs_close[trainingProbes,],coefEsts)


## close et al single cell
pdf("closeSingle_compEsts.pdf", w=50,h=6)
palette(pal)
par(mar = c(9,6,2,2), cex.axis=1.5,cex.lab=2)
ooClose= order(rse_gene_close$SampleLabel)
bp=barplot(t(close_PropEsts[ooClose,]), col = pal,
	ylim = c(0,1.3),xaxt = "n", yaxt= "n",ylab="Class Proportion")

g = split(bp,rse_gene_close$SampleLabel[ooClose])
axis(1,at=sapply(g,mean), font=2, labels = gsub(" day ", "", names(g)),las=3)
abline(v = sapply(g,min)[-1]-0.5,lwd=4,col="black")
legend("top", colnames(close_PropEsts), pch = 15, 
	col = 1:8,bg="white", nc = 4,cex=1.5,pt.cex=2)
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
rse_geneCamp$SampleLabel = paste0(rse_geneCamp$Stage_s, "_", rse_geneCamp$Day)
yExprs_camp = log2(getRPKM(rse_geneCamp,"Length")+1)
camp_PropEsts= minfi:::projectCellType(yExprs_camp[trainingProbes,],coefEsts)


## camp et al data
pdf("campFetal_compEsts.pdf", w=50,h=6)
palette(pal)
par(mar = c(9,6,2,2), cex.axis=1.5,cex.lab=2)
ooCamp= order(rse_geneCamp$source_name_s,rse_geneCamp$Day)
bp=barplot(t(camp_PropEsts[ooCamp,]), col = pal,
	ylim = c(0,1.3),xaxt = "n", yaxt= "n",ylab="Class Proportion")

g = split(bp,rse_geneCamp$SampleLabel[ooCamp])
axis(1,at=sapply(g,mean), font=2, labels = gsub(" day ", "", names(g)),las=3)
abline(v = sapply(g,min)[-1]-0.5,lwd=4,col="black")
legend("top", colnames(camp_PropEsts), pch = 15, 
	col = 1:8,bg="white", nc = 4,cex=1.5,pt.cex=2)
axis(2, at=seq(0,1,by=0.25))
dev.off()

#####################
## organoids ########
library(GEOquery)

load("/dcl01/ajaffe/data/lab/singleCell/pasca_organoid/rse_gene_Pasca_Organoids_human_n730.Rdata")
rse_genePasca = rse_gene
yExprs_pasca = log2(getRPKM(rse_genePasca,"Length")+1)
pasca_PropEsts= minfi:::projectCellType(yExprs_pasca[rownames(coefEsts),],coefEsts)

pdPasca = read.delim("/dcl01/ajaffe/data/lab/singleCell/pasca_organoid/pasca_SraRunTable.txt",
	as.is=TRUE)
pdPasca = pdPasca[match(colnames(rse_genePasca), pdPasca$Run),]

rse_genePasca$Age = as.numeric(ss(pdPasca$age, " ", 2))
rse_genePasca$CellType = factor(pdPasca$cell_type)
rse_genePasca$SampleLabel = paste0(rse_genePasca$CellType , "_", rse_genePasca$Age)

## camp et al data
pdf("pascaOrganoid_compEsts.pdf", w=50,h=6)
palette(pal)
par(mar = c(9,6,2,2), cex.axis=1.5,cex.lab=2)
ooPasca = order(rse_genePasca$CellType,rse_genePasca$Age)
bp=barplot(t(pasca_PropEsts[ooPasca,]), col = pal,
	ylim = c(0,1.3),xaxt = "n", yaxt= "n",ylab="Class Proportion")

g = split(bp,rse_genePasca$SampleLabel[ooPasca])
axis(1,at=sapply(g,mean), font=2, labels = gsub(" day ", "", names(g)),las=3)
abline(v = sapply(g,min)[-1]-0.5,lwd=4,col="black")
legend("top", colnames(pasca_PropEsts), pch = 15, 
	col = 1:8,bg="white", nc = 4,cex=1.5,pt.cex=2)
axis(2, at=seq(0,1,by=0.25))
dev.off()