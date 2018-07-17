#####

library(SummarizedExperiment)
library(recount)
library(jaffelab)
library(RColorBrewer)
library(minfi)

## libd time course
MAINDIR = "/dcl01/lieber/ajaffe/lab/libd_stem_timecourse"
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

## human data from span
load("/dcl01/lieber/ajaffe/lab/brainspan_analysis/rawCounts_full_brainspan_n607.rda", verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/brainspan_analysis/full.brainspan_final.phenotype.rda",verbose=TRUE)
pdSpan = cbind(fin_pdSpan, metrics[fin_pdSpan$lab,])
pdSpan$wig = NULL
geneCounts = geneCounts[,pdSpan$lab]
rownames(pdSpan) = pdSpan$lab
rm(exonCounts, exonMap, jCounts, jMap)

## Create gene,exon RangedSummarizedExperiment objects
gr_genes <- GRanges(seqnames = geneMap$Chr,
    IRanges(geneMap$Start, geneMap$End), strand = geneMap$Strand)
names(gr_genes) <- rownames(geneMap)
mcols(gr_genes) <- DataFrame(geneMap[, - which(colnames(geneMap) %in%
    c('Chr', 'Start', 'End', 'Strand'))])
rse_geneSpan <- SummarizedExperiment(assays = list('counts' = geneCounts),
    rowRanges = gr_genes, colData = pdSpan)

###################
## based on stem ##
###################
geneRpkmTC = getRPKM(rse_geneTC, "Length")
yRpkmTC = log2(geneRpkmTC+1)
pcaTC = prcomp(t(yRpkmTC))
pcaTC_Vars = getPcaVars(pcaTC)

# project
geneRpkmSpan = getRPKM(rse_geneSpan,length_var="Length")
yExprsSpan_Scaled = scale(t(log2(geneRpkmSpan+1)), pcaTC$center, pcaTC$scale) 
genePCs_Span_projected = yExprsSpan_Scaled %*% pcaTC$rotation 


## panel A
pdf("brainspan_projection_PCsbasedOnLIBD.pdf", useDingbats=FALSE, h=6,w=7)
par(mar=c(5,6,4,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(9,"Spectral"))
plot(pcaTC$x, bg = as.numeric(factor(rse_geneTC$COND))+3,pch=21,cex=2,
	xlab=paste0("PC1: ", pcaTC_Vars[1], "% Var Expl"),
    ylab=paste0("PC2: ", pcaTC_Vars[2], "% Var Expl"))
legend("bottom", paste0(levels(factor(rse_geneTC$COND))), col=4:9,
       pch = 16, cex=.85)
	   
palette(c(colorRampPalette(c("red","white"))(11), 
	colorRampPalette(c("white", "blue"))(19)))
points(genePCs_Span_projected[,1:2], pch = 22, cex=2,
	bg=as.numeric(factor(rse_geneSpan$Age)))
legend("topright", c("Age: -0.615","Birth","Age: 40"), col=c(1,"black",30), 
       pch = c(15,22,15), cex=.9)
dev.off()
	   

#########################
# based on brain  ##
#########################

yRpkmSpan= log2(geneRpkmSpan+1)
pcaSpan = prcomp(t(yRpkmSpan))
pcaSpan_Vars = getPcaVars(pcaSpan)

# project
yExprsTC_Scaled = scale(t(yRpkmTC), pcaSpan$center, pcaSpan$scale) 
genePCs_TC_projected = yExprsTC_Scaled %*% pcaSpan$rotation 

## panel A
pdf("brainspan_projection_PCsbasedOnBrainSpan.pdf", 
	useDingbats=FALSE, h=6,w=7)
par(mar=c(5,6,4,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(c(colorRampPalette(c("red","white"))(11), 
	colorRampPalette(c("white", "blue"))(19)))
plot(pcaSpan$x, bg = as.numeric(factor(rse_geneSpan$Age)), 
	cex=2, pch=22, 
	xlab=paste0("PC1: ", pcaSpan_Vars[1], "% Var Expl"),
    ylab=paste0("PC2: ", pcaSpan_Vars[2], "% Var Expl"))
legend("bottom", c("Age: -0.615","Birth","Age: 40"), col=c(1,"black",30), 
       pch = c(15,22,15), cex=.8)
palette(brewer.pal(9,"Spectral"))
points(genePCs_TC_projected[,1:2], pch = 21, cex=2,
	bg=as.numeric(factor(rse_geneTC$COND))+3)
legend("bottomleft", paste0(levels(factor(rse_geneTC$COND))), col=4:9,
       pch = 16, cex=.8)
dev.off()



## pc1 alone??
# rse_geneTC$CONDITION = factor(rse_geneTC$CONDITION, 
	# levels = c("RENEW", "ACC_DORSAL(2)", "NPC", "ROSETTE", 
		# "NEURONS_ALONE", "NEURONS_PLUS_ASTROS"))
par(mar=c(9,6,2,2))
boxplot(genePCs_TC_projected[,1] ~ rse_geneTC$COND,
	xlim = c(0.5,38.5), ylim = c(-100,100),las=3,
	ylab= paste0("PC1: ", pcaSpan_Vars[1], "% Var Expl"))
boxplot(pcaSpan$x[,1] ~ factor(round(rse_geneSpan$Age,2)), 
	add=TRUE, at = 8:37, las = 3)	
abline(v=7,lty=2)



## pc1 alone??
pdf("brainspan_boxplot_projection_PCsbasedOnBrainSpan.pdf", useDingbats=FALSE, h=6,w=25)
# boxplot(genePCs_TC_projected[,1] ~ as.numeric(ss(rse_geneTC$DAY,"_",2)))
# rse_geneTC$COND = factor(rse_geneTC$COND, 
	# levels = c("RENEW", "ACC_DORSAL(2)", "NPC", "ROSETTE", 
		# "NEURONS_ALONE", "NEURONS_PLUS_ASTROS"))
par(mar=c(10,6,2,2),cex.axis=1.5,cex.lab=2,cex.main=2)
palette(brewer.pal(9,"Spectral"))
boxplot(-genePCs_TC_projected[,1] ~ rse_geneTC$COND, outline=FALSE, xaxt="n",
		xlim = c(1.5,36.5), ylim = c(-110,140), las=3,
		ylab= paste0("PC1: ", pcaSpan_Vars[1], "% Var Expl"))
	points(-genePCs_TC_projected[,1] ~ jitter(as.numeric(rse_geneTC$COND),amount=0.15), pch = 21, 
		bg=as.numeric(factor(rse_geneTC$COND))+3, cex=2)
	axis(1, at=1:6, labels = FALSE)
	text(1.1:6.1, par("usr")[3]-12,
		srt=40, adj= 1, xpd = TRUE,
		labels = levels(rse_geneTC$COND), cex=1.25)
palette(c(colorRampPalette(c("red","white"))(11), 
	colorRampPalette(c("white", "blue"))(19)))
boxplot(-pcaSpan$x[,1] ~ round(rse_geneSpan$Age,2), outline=FALSE, yaxt="n",
	at = 8:37, add=TRUE, las=3)	
	points(-pcaSpan$x[,1] ~ jitter(as.numeric(factor(rse_geneSpan$Age))+7,amount=0.15), pch = 22, 
		bg=as.numeric(factor(rse_geneSpan$Age)), cex=2)

abline(v=7,lty=2)

dev.off()



### using AgeCat
rse_geneSpan$AgeCat = "Adult"
rse_geneSpan$AgeCat[rse_geneSpan$Age<20] = "Teens"
rse_geneSpan$AgeCat[rse_geneSpan$Age<10] = "Child"
rse_geneSpan$AgeCat[rse_geneSpan$Age<1] = "Infant"
rse_geneSpan$AgeCat[rse_geneSpan$Age<0] = "Late Fetal"
rse_geneSpan$AgeCat[rse_geneSpan$Age<(-.25)] = "Mid Fetal"
rse_geneSpan$AgeCat[rse_geneSpan$Age<(-.5)] = "Early Fetal"
rse_geneSpan$AgeCat = as.factor(rse_geneSpan$AgeCat)
rse_geneSpan$AgeCat = factor(rse_geneSpan$AgeCat,levels(rse_geneSpan$AgeCat)[c(3,6,5,4,2,7,1)])


pdf("brainspan_boxplot_projection_PCsbasedOnBrainSpan_AgeCat.pdf", useDingbats=FALSE, h=6,w=10)
par(mar=c(8,6,2,2),cex.axis=1.5,cex.lab=2,cex.main=2)
palette(brewer.pal(9,"Spectral"))
boxplot(-genePCs_TC_projected[,1] ~ rse_geneTC$COND, outline=FALSE, xaxt="n",
		xlim = c(0.5,14.5), ylim = c(-110,140), las=3,
		ylab= paste0("PC1: ", pcaSpan_Vars[1], "% Var Expl"))
	points(-genePCs_TC_projected[,1] ~ jitter(as.numeric(rse_geneTC$COND),amount=0.15), pch = 21, 
		bg=as.numeric(factor(rse_geneTC$COND))+3, cex=2)
	axis(1, at=1:6, labels = FALSE)
	text(1.1:6.1, par("usr")[3]-12,
		srt=40, adj= 1, xpd = TRUE,
		labels = levels(rse_geneTC$COND), cex=1.25)
# palette(brewer.pal(9,"RdBu"))
palette(c(colorRampPalette(c("red","white"))(11), 
	colorRampPalette(c("white", "blue"))(19)))
boxplot(-pcaSpan$x[,1] ~ rse_geneSpan$AgeCat, outline=FALSE, yaxt="n", xaxt="n",
	at = 8:14, add=TRUE, las=3)	
	points(-pcaSpan$x[,1] ~ jitter(as.numeric(factor(rse_geneSpan$AgeCat))+7,amount=0.15), pch = 22, 
		bg=as.numeric(factor(rse_geneSpan$AgeCat))*4-2, cex=2)
	axis(1, at=8:14, labels = FALSE)
	text(8.1:14.1, par("usr")[3]-12,
		srt=40, adj= 1, xpd = TRUE,
		labels = levels(rse_geneSpan$AgeCat), cex=1.25)
abline(v=7,lty=2)

dev.off()

##############################################
###### regression calibration approach  ######

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

## brainspan NCX
rse_geneSpan$RegionGroup = rse_geneSpan$Regioncode 
rse_geneSpan$RegionGroup[rse_geneSpan$Regioncode %in% 
				c("A1C", "DFC", "IPC","ITC", "M1C", "MFC",
				"OFC", "S1C", "STC", "V1C", "VFC")] = "NCX"
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
numProbes=20
probeList <- lapply(tstatList, function(x) {
	y <- x[x[, "p.value"] < 1e-20, ]
	yUp <- y[order(y[, "dm"], decreasing = TRUE), ]
	yDown <- y[order(y[, "dm"], decreasing = FALSE),
		]
	c(rownames(yUp)[1:numProbes], rownames(yDown)[1:numProbes])
})

## filter
trainingProbes <- unique(unlist(probeList))
mergeMarkerExprs <- yExprs_Merge[trainingProbes, ]
mergeMarkerMeanExprs <- colMeans(mergeMarkerExprs)
names(mergeMarkerMeanExprs) <- names(tIndexes)

form <- as.formula(sprintf("y ~ %s - 1", paste(levels(group),collapse = "+")))
phenoDF <- as.data.frame(model.matrix(~group - 1))
colnames(phenoDF) <- sub("^group", "", colnames(phenoDF))

## do calibration
coefEsts <- minfi:::validationCellType(Y = mergeMarkerExprs, 
	pheno = phenoDF, modelFix = form)$coefEsts
save(coefEsts, file = "iPSC_BrainSpan_coefEsts_calibration.rda")
load("iPSC_BrainSpan_coefEsts_calibration.rda")
write.csv(coefEsts, file="iPSC_BrainSpan_coefEsts_calibration.csv")

## save output as text
signif(coefEsts[probeList$NCX_Adult,],3)

#######################
### fit model
stemPropEsts = minfi:::projectCellType(yRpkmTC[rownames(coefEsts),],
	coefEsts,lessThanOne=TRUE)

## test for differences
modAdd = model.matrix(~ as.numeric(COND) + LINE + totalAssignedGene, 
	data=colData(rse_geneTC))
t(apply(stemPropEsts, 2, function(y) summary(lm(y ~ modAdd - 1))$coef[2,]))
# > t(apply(stemPropEsts, 2, function(y) summary(lm(y ~ mod - 1))$coef[2,]))
                    # Estimate   Std. Error     t value     Pr(>|t|)
# iPSC           -1.222163e-01 7.202905e-03 -16.9676476 7.450692e-33
# NCX_EarlyFetal  4.021083e-02 6.695292e-03   6.0058366 2.364852e-08
# NCX_MidFetal    4.618401e-02 3.905727e-03  11.8246897 1.737101e-21
# NCX_LateFetal   1.367008e-02 2.797509e-03   4.8865203 3.414477e-06
# NCX_Infant      8.095303e-03 1.445778e-03   5.5992699 1.530725e-07
# NCX_Child      -3.335493e-17 3.890172e-18  -8.5741534 5.948944e-14
# NCX_Teens      -6.052086e-18 6.132057e-18  -0.9869585 3.257716e-01
# NCX_Adult       2.169578e-02 1.666090e-03  13.0219704 3.111767e-24
mod = model.matrix(~ COND + LINE + totalAssignedGene, 
	data=colData(rse_geneTC))
statList = lapply(as.data.frame(stemPropEsts), function(y) summary(lm(y ~ mod - 1))$coef[2:6,])
tMat = t(sapply(statList, function(x) x[,3]))
pMat = t(sapply(statList, function(x) x[,4]))
colnames(tMat) = colnames(pMat) = gsub("modCOND", "", colnames(tMat))

round(tMat,2)
signif(pMat,3)

## mean proportions
cIndexes =splitit(rse_geneTC$COND)
meanMat = round(100*t(sapply(cIndexes, function(ii) colMeans(stemPropEsts[ii,]))),2)
library(pheatmap)
pal = colorRampPalette(c("white", brewer.pal(n = 7, name = "Reds")))(100)
pheatmap(meanMat[-1,c(1:5,8)], color = pal, cluster_rows =FALSE, cluster_cols = FALSE)

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

##############################
### validate using BrainSeq
load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseGene_merged_n732.rda")
rse_geneDLPFC = rse_gene[,rse_gene$Dx == "Control" & 
	rse_gene$Race %in% c("AA","CAUC")]
rse_geneDLPFC$ageGroup = cut(rse_geneDLPFC$Age, c(-1,0,1,10,20,50,100))
levels(rse_geneDLPFC$ageGroup) = c("Fetal","Infant","Child","Teens","Adult","50+")

rse_geneDLPFC = rse_geneDLPFC[,order(rse_geneDLPFC$Age)]
yExprsDLPFC = log2(getRPKM(rse_geneDLPFC,length_var="Length")+1)

# project
dlpfcPropEsts = minfi:::projectCellType(yExprsDLPFC[trainingProbes,],coefEsts)

## visualize
library(pheatmap)
 col.pal = brewer.pal(9,"Blues")

pdf("comp_heatmap.pdf", useDingbats=FALSE, w=8,h=16)
rownames(stemPropEsts) = rse_geneTC$COND
pheatmap(stemPropEsts, cluster_rows=TRUE, 
		cluster_cols=FALSE, color=col.pal)
rownames(dlpfcPropEsts) = rse_geneDLPFC$ageGroup
pheatmap(dlpfcPropEsts, cluster_rows=FALSE, 
		cluster_cols=FALSE, color=col.pal)
dev.off()

#### comp bar charts
pal = c("#762a83", brewer.pal(7,"RdYlBu"))


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