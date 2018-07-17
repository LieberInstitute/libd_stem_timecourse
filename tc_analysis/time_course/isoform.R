### load packages
library(jaffelab)
library(GenomicRanges)
library(rtracklayer)
library(SummarizedExperiment)
library(RColorBrewer)
library(rafalib)
library(edgeR)
library(limma)
library(recount)


#############################
#### load feature counts ####
MAINDIR = '/dcl01/lieber/ajaffe/lab/libd_stem_timecourse'
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n146.rda"))
load(file.path(MAINDIR,"data/libd_stemcell_rseTx_counts.rda"))


########################
#### subset samples ####
dropInd = which(rse_gene$Class != "Naked genomes" | rse_gene$CONDITION %in% c("RENEW","NEURONS_ALONE"))

rse_gene = rse_gene[,-dropInd]  ## n=106
rse_tx = rse_tx[,-dropInd]

## Recommended cutoffs from Leo's function
# Gene
# 2017-12-13 13:52:30 the suggested expression cutoff is 0.25 # 0.25 0.21
# Exon
# 2017-12-13 13:58:29 the suggested expression cutoff is 0.27 # 0.27 0.22
# Jxn
# 2017-12-13 14:01:13 the suggested expression cutoff is 0.36 # 0.29 0.36
# Tx
# 2017-12-13 14:42:10 the suggested expression cutoff is 0.31 # 0.31 0.24

###########################
#### filter expression ####
geneRpkm = recount::getRPKM(rse_gene, length_var="Length")
txTpm = assays(rse_tx)$tpm

mcols(rse_gene)$meanExprs = rowMeans(geneRpkm)
rse_gene = rse_gene[which(mcols(rse_gene)$meanExprs > 0.1),]

# mcols(rse_tx)$meanExprs = rowMeans(txTpm)
# rse_tx = rse_tx[which(mcols(rse_tx)$meanExprs > 0.3),]


########################
#### counts objects ####
geneCounts = assays(rse_gene)$counts
geneMap = as.data.frame(rowRanges(rse_gene))
txCounts = assays(rse_tx)$counts
txMap = as.data.frame(rowData(rse_tx))

##################
#### clean pd ####
pd = colData(rse_gene)
pd$REP = as.numeric(pd$BIO_REP)
pd$REP[is.na(pd$REP)] = 1
pd$DAY = ordered(pd$DAY, levels=sort(as.numeric(levels(as.factor(pd$DAY)))) )
pd$DIV = as.numeric(as.character(pd$DAY))
pd$EXPERIMENT[pd$EXPERIMENT=="TIMECOURSE1"] = "TIMECOURSE"
pd$LINE = gsub('-','_',pd$LINE)

pd$COND = pd$CONDITION
pd$COND[grep("NEURON",pd$COND)] = "NEURON"  ## rename NEURONS_PLUS_ASTROS
pd$COND[grep("ACC_DORSAL",pd$COND)] = "ACC_DORSAL"  ## rename ACC_DORSAL(2)
pd$COND = as.factor(pd$COND)
pd$COND = factor(pd$COND,levels(pd$COND)[c(1,3,4,2)])  ## put levels in order


load("voomStats_4features.rda", verbose=TRUE)


###############
#### plots ####
###############

# Colors by Donor
pal = c(colorRampPalette(brewer.pal(10,"Paired")[1:2])(4), 
		colorRampPalette(brewer.pal(10,"Paired")[3:4])(3), 
		brewer.pal(10,"Paired")[9:10], brewer.pal(10,"Paired")[5:8])
pd$lineCol = pal[as.numeric(factor(pd$LINE))]

geneRpkm = recount::getRPKM(rse_gene, length_var="Length")
txTpm = assays(rse_tx)$tpm
gExprs = as.matrix(log2(geneRpkm+1))
txExprs = as.matrix(log2(txTpm+1))


#######################
#### major isoform ####

## average txTpm for each condition
txAvg = data.frame(matrix(NA, nrow = nrow(txTpm), ncol = 4))
rownames(txAvg) = rownames(txTpm)
colnames(txAvg) = levels(pd$COND)
for (i in 1:4) {
	txAvg[,i] = rowMeans(txTpm[,which(pd$COND==levels(pd$COND)[i])])
}

## top tx of for each gene at each condition (25466 x 4)
majorIsoform = data.frame(matrix(NA, nrow = nrow(geneRpkm), ncol = 4))		# txID
majorIsoPercent = data.frame(matrix(NA, nrow = nrow(geneRpkm), ncol = 4))	# percent
majorIsoTPM = data.frame(matrix(NA, nrow = nrow(geneRpkm), ncol = 4))		# TPM
majorIsoAmountP = data.frame(matrix(NA, nrow = nrow(geneRpkm), ncol = 4))	# Diff in percent
majorIsoAmountTPM = data.frame(matrix(NA, nrow = nrow(geneRpkm), ncol = 4))	# Diff in TPM

rownames(majorIsoform) = rownames(majorIsoPercent) = rownames(majorIsoTPM) = rownames(majorIsoAmountP) = rownames(majorIsoAmountTPM) = rownames(geneRpkm)
colnames(majorIsoform) = colnames(majorIsoPercent) = colnames(majorIsoTPM) = colnames(majorIsoAmountP) = colnames(majorIsoAmountTPM) = levels(pd$COND)


## list of data.frames with tx info for each gene
txTpms = list()
txPercents = list()

## load below
# # for (i in 1:nrow(geneMap)) {
	# # gene = rownames(geneMap)[i]
	# # txs = geneMap$gencodeTx[[i]]
	# # tmp = txAvg[txs,]

	# # txTpms[[i]] = tmp
	# # txPercents[[i]] = t(t(tmp)/rowSums(t(tmp)))
	# # names(txTpms)[i] = names(txPercents)[i] = gene
	
	# # majorIsoPercent[i,] = apply(txPercents[[i]], 2, max)
	# # majorIsoTPM[i,] = apply(txTpms[[i]], 2, max)
	
	# # if (nrow(tmp) == 1) {	
		# # majorIsoAmountP[i,] = 1
		# # majorIsoAmountTPM[i,] = 1
	# # } else {
		# # ## 2nd highest percentage
		# # majorIsoAmountP[i,] = apply(txPercents[[i]], 2, function(x){sort(x)[nrow(txPercents[[i]])-1][1]})
		# # ## diff between major iso and next highest percentage
		# # majorIsoAmountP[i,] = majorIsoPercent[i,] - majorIsoAmountP[i,]
		
		# # ## 2nd highest tpm
		# # majorIsoAmountTPM[i,] = apply(txTpms[[i]], 2, function(x){sort(x)[nrow(txTpms[[i]])-1][1]})
		# # ## diff between major iso and next highest tpm
		# # majorIsoAmountTPM[i,] = majorIsoTPM[i,] - majorIsoAmountTPM[i,]
	# # }

	# # ## top tx at each condition
	# # for (j in 1:4) {
		# # majorIsoform[i,j] = ifelse(max(tmp[,j])==0,NA,rownames(tmp)[which(tmp[,j] == max(tmp[,j]))])
	# # }
# # }

# save(majorIsoform,majorIsoPercent,majorIsoTPM,majorIsoAmountP,majorIsoAmountTPM,
		# txTpms,txPercents, file="isoforms.RData")

load("isoforms.RData")

shift = rep("",nrow(majorIsoform) )
for (i in 1:nrow(majorIsoform)) {
	## not a shift if only 1 unique tx (among non-zero conditions)
	shift[i] = ifelse(length(unique(as.character(majorIsoform[i,!is.na(majorIsoform[i,1:4])])))==1,"","shift")
	shift[i] = ifelse(nrow(txTpms[[i]])==1,"",shift[i]) ## no shift if only 1 tx
}
shiftInd = which(shift == "shift")

table(shift)
# 		 shift
# 17586  7880


majorIsoform2 = majorIsoform[shiftInd,]
majorIsoAmountTPM2 = majorIsoAmountTPM[shiftInd,]
txTpms2 = txTpms[shiftInd]
majorIsoAmountP2 = majorIsoAmountP[shiftInd,]
majorIsoPercent2 = majorIsoPercent[shiftInd,]
# txPercents2 = txPercents[shiftInd]

## shifts in DE genes
sum(rownames(majorIsoform2) %in% rownames(gVoomStats)[which(gVoomStats$adj.P.Val < 0.01)])
# [1] 7108

## shifts between AD and NPC
length(which(majorIsoform2$ACC_DORSAL != majorIsoform2$NPC))
# [1] 2390
## DE genes
sum(rownames(majorIsoform2)[which(majorIsoform2$ACC_DORSAL != majorIsoform2$NPC)] %in% rownames(gVoomStats)[which(gVoomStats$"q_NPC-ACC_DORSAL" < 0.01)])
# [1] 1216

## shifts between NPC and ROSE
length(which(majorIsoform2$ROSETTE != majorIsoform2$NPC))
# [1] 1989
## DE genes
sum(rownames(majorIsoform2)[which(majorIsoform2$ROSETTE != majorIsoform2$NPC)] %in% rownames(gVoomStats)[which(gVoomStats$"q_ROSETTE-NPC" < 0.01)])
# [1] 238

## shifts between ROSE and NEURON
length(which(majorIsoform2$ROSETTE != majorIsoform2$NEURON))
# [1] 6479
## DE genes
sum(rownames(majorIsoform2)[which(majorIsoform2$ROSETTE != majorIsoform2$NEURON)] %in% rownames(gVoomStats)[which(gVoomStats$"q_NEURON-ROSETTE" < 0.01)])
# [1] 4011





##############################
######## load in txdb ########
### (prepare for plotting) ###
library(GenomicFeatures)
## GENCODE
chrInfo = getChromInfoFromUCSC("hg38")
chrInfo$chrom = as.character(chrInfo$chrom)
chrInfo = chrInfo[chrInfo$chrom %in% paste0("chr", c(1:22,"X","Y", "M")),]
chrInfo$isCircular = rep(c(FALSE, TRUE), c(24,1))
si = with(chrInfo, Seqinfo(as.character(chrom), length, isCircular, genome="hg38"))

gencode_v25 = import(con = "/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf", format = "gtf")
seqlevels(gencode_v25, pruning.mode="coarse") = paste0("chr", c(1:22,"X","Y"))
seqinfo(gencode_v25) = si
gencode_v25_txdb = makeTxDbFromGRanges(gencode_v25)
tx = exonsBy(gencode_v25_txdb, use.names=TRUE)


# gencode_v25_37 = import(con = "/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh37_hg19/gencode.v25lift37.annotation.gtf")
# seqlevels(gencode_v25_37, pruning.mode="coarse") = "chr10"
# gencode_v25_37_txdb = makeTxDbFromGRanges(gencode_v25_37)
# select(gencode_v25_37_txdb, keys="ENSG00000177283.6_2", columns="TXNAME", keytype="GENEID")
# 'select()' returned 1:many mapping between keys and columns
               # GENEID              TXNAME
# 1 ENSG00000177283.6_2 ENST00000374694.2_1
# 2 ENSG00000177283.6_2 ENST00000579659.1_1

# remove seqlengths so the extra intron doesn't get plotted
tmp = seqinfo(tx)
tmp@seqlengths <- as.integer(rep(NA, 25))
seqinfo(tx) = tmp
##############################

## edit derfinderPlot:::.plotPolygon for diff colors per tx
.plotPolygonColorsExon = function (info) {
	exonsPerTX = unlist(lapply(tx[subjectHits(ov)][keepInd], length))
	exonsInd = data.frame(exonsPerTX, endInd = cumsum(exonsPerTX))
	exonsInd$startInd = exonsInd$endInd - exonsInd$exonsPerTX + 1
	for (trx in 1:nrow(exonsInd)) {
		for (row in exonsInd$startInd[trx]:exonsInd$endInd[trx]) {
			polygon(info$x[row, ], info$y[row, ], lwd = 1.2, col = gg_colors[trx])
		}
	}
}


#############
## tx lengths
txLengths = as.data.frame(gencode_v25[which(mcols(gencode_v25)$type=="transcript")])[,c(1:7,10,13,16)]
rownames(txLengths) = txLengths$transcript_id
txLengths = split(txLengths, txLengths$gene_id)  ## turn into list by gene name
txLengths = txLengths[names(txTpms)]  ## put in order of txTpms
txLengths = txLengths[rownames(majorIsoform2)]

## check number of switches that have a tx < 250 bp
txLengthsIsoform = list()
for (i in 1:nrow(majorIsoform2) ) {
	uniquetx = unique(unlist(majorIsoform2[i,]))
	txLengthsIsoform[[i]] = txLengths[[i]][uniquetx,]
}
names(txLengthsIsoform) = rownames(majorIsoform2)
txSmall = lapply(txLengthsIsoform, function(x) min(x$width) )
table(txSmall < 250)
# FALSE  TRUE
#  7841    24

# ## get proportion expressed on log scale
# txTpmLog = lapply(txTpms2, function(x) log2(x+1) )
# percentLog = lapply(txTpmLog, function(x)  t(t(x)/rowSums(t(x))) )
# percentLogIso = t(as.data.frame(lapply(percentLog, function(x) apply(x, 2, max) )))

############################################
#####    plot tx tpms and structure   ######
library(derfinderPlot)
library(GenomicFeatures)

gg_colors = c("#a38692","#70b670","#3290d8","#be9f52","#b12a2d","#7c52be","#df9c7f","#5e756b","#bbbed2",
		"#b8c85e","#7a633e","#fb895a","#4d5382")

## number to plot:
n = 100
##
majorIsoAmountP2 = as.data.frame(majorIsoAmountP2[order(rowMeans(majorIsoAmountP2),decreasing=TRUE),])
majorIsoAmountP2$smallTxLength = txSmall[rownames(majorIsoAmountP2)]
# majorIsoAmountP2 = majorIsoAmountP2[-which(majorIsoAmountP2$smallTxLength <250) , ]  ## uncomment to drop <250bp
top6 = rownames(majorIsoAmountP2)[1:n]
top6tpm = txTpms2[top6]
top6stats = gVoomStats[top6,]

pdf("isoform_top_shifts_with_structure.pdf", h=10, w=12, onefile=TRUE)
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
for (i in 1:n) {
	genei = paste0(rownames(top6stats)[i],"\n",top6stats$Symbol[i])
	tpmi = top6tpm[[i]]
	nTx = nrow(tpmi)

	tpmMat = as.matrix(log2(tpmi+1))
	perMat = t(t(tpmMat)/rowSums(t(tpmMat)))
	colnames(tpmMat) = colnames(perMat) = c("Acc Dorsal","NPC","Rosette","Neuron")
	tpmMat[which(tpmMat == 0, arr.ind=TRUE)] = 0.01
	perMat[which(perMat == 0, arr.ind=TRUE)] = 0.01

	par(mar=c(3,6,2,2))
	## first barplot: log(tpm)
	barplot(tpmMat, beside=T, col=gg_colors[1:nrow(tpmMat)],
			cex.axis=1.4, ylim=c(0,max(tpmMat)*1.05), cex.names=1.5, cex.lab=1.5,
			ylab = "log2(TPM+1)")
	box()
	## second barplot: proportions
	barplot(perMat, beside=T, col=gg_colors[1:nrow(perMat)], 
			cex.axis=1.4, ylim=c(0,1.01), cex.names=1.5, cex.lab=1.5,
			ylab = "Relative proportion expressed")
	box()

	par(mar=c(6,15,5,10))
	## plot tx exon structure
	gene = top6[i]
	geneInfo = geneMap[gene,]
	txsOfGene = unlist(geneInfo$gencodeTx)
	chrom = geneInfo$seqnames
	minRange = geneInfo$start
	maxRange = geneInfo$end
	strnd=geneInfo$strand
	gr = GRanges(chrom, IRanges(minRange,maxRange), strnd)
	yLabelCoords = -1:-length(txsOfGene)  ## negative y-coords so first tx is on top

	## plot
	ov <- findOverlaps(gr, tx)
	## drop txs in range that aren't part of gene
	txInOv = names(tx[subjectHits(ov)]) ##
	keepInd = which(txInOv %in% txsOfGene) ##
	txList <- split(tx[subjectHits(ov)][keepInd], queryHits(ov)[keepInd])
	txIDs <- names(tx[subjectHits(ov)])[keepInd]
	poly.data <- lapply(txList, derfinderPlot:::.plotData)[[1]]
	if(strnd=="+") { 
		poly.data$intron$y = - poly.data$intron$y
		poly.data$exon$y = - poly.data$exon$y	}  ## negative so first is on top 
	yrange <- range(poly.data$exon$y)
	xrange <- c(start(gr),end(gr))
	plot(0,0, type="n", main=genei, xlim = xrange , ylim = yrange + c(-0.75, 0.75),
		yaxt="n", ylab="", xlab=as.character(seqnames(gr)),
		cex.axis = 1.3, cex.lab = 1.3, cex.main = 2)
	axis(2, at=yLabelCoords, labels=txIDs, las=1, cex.axis = 1.5)
	derfinderPlot:::.plotPolygon(poly.data$intron, 'gray90')
	.plotPolygonColorsExon(poly.data$exon)

}
dev.off()





###### ###### ###### ###### 
###### New numbers after taking out 24 short transcripts

geneSmall = names(txSmall)[which(txSmall<250)]
smallInd = which(rownames(majorIsoform2) %in% geneSmall)
majorIsoform3 = majorIsoform2[-smallInd,]

## shifts in DE genes
sum(rownames(majorIsoform3) %in% rownames(gVoomStats)[which(gVoomStats$adj.P.Val < 0.01)])
# [1] 7086

## shifts between AD and NPC
length(which(majorIsoform3$ACC_DORSAL != majorIsoform3$NPC))
# [1] 2383
## DE genes
sum(rownames(majorIsoform3)[which(majorIsoform3$ACC_DORSAL != majorIsoform3$NPC)] %in% rownames(gVoomStats)[which(gVoomStats$"q_NPC-ACC_DORSAL" < 0.01)])
# [1] 1213

## shifts between NPC and ROSE
length(which(majorIsoform3$ROSETTE != majorIsoform3$NPC))
# [1] 1982
## DE genes
sum(rownames(majorIsoform3)[which(majorIsoform3$ROSETTE != majorIsoform3$NPC)] %in% rownames(gVoomStats)[which(gVoomStats$"q_ROSETTE-NPC" < 0.01)])
# [1] 238

## shifts between ROSE and NEURON
length(which(majorIsoform3$ROSETTE != majorIsoform3$NEURON))
# [1] 6459
## DE genes
sum(rownames(majorIsoform3)[which(majorIsoform3$ROSETTE != majorIsoform3$NEURON)] %in% rownames(gVoomStats)[which(gVoomStats$"q_NEURON-ROSETTE" < 0.01)])
# [1] 3998


###### ###### ###### ###### 
###### Summary of results
# 7856 out of 25,466 genes have a shift in major isoform
# 7086 out of 20,220 DE genes (FDR < 0.01) have a shift

# 2383 shift between AD and NPC
# 1213 shift between 9067 DE genes AD and NPC

# 1982 shift between NPC and ROSETTE
# 238 shift between 1994 DE genes NPC and ROSETTE

# 6459 shift between ROSETTE and NEURON
# 3998 shift between 12,951 DE genes ROSETTE and NEURON
###### ###### ###### ###### 


























##############################
##### plot tx structure ######
library(derfinderPlot)
library(GenomicFeatures)

######## load in txdb ########
## GENCODE
chrInfo = getChromInfoFromUCSC("hg38")
chrInfo$chrom = as.character(chrInfo$chrom)
chrInfo = chrInfo[chrInfo$chrom %in% paste0("chr", c(1:22,"X","Y", "M")),]
chrInfo$isCircular = rep(c(FALSE, TRUE), c(24,1))
si = with(chrInfo, Seqinfo(as.character(chrom), length, isCircular, genome="hg38"))

gencode_v25 = import(con = "/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf", format = "gtf")
seqlevels(gencode_v25, pruning.mode="coarse") = paste0("chr", c(1:22,"X","Y","M"))
seqinfo(gencode_v25) = si
gencode_v25_txdb = makeTxDbFromGRanges(gencode_v25)
tx = exonsBy(gencode_v25_txdb, use.names=TRUE)

# remove seqlengths so the extra intron doesn't get plotted
tmp = seqinfo(tx)
tmp@seqlengths <- as.integer(rep(NA, 25))
seqinfo(tx) = tmp
##############################

## edit derfinderPlot:::.plotPolygon for diff colors per tx
.plotPolygonColorsExon = function (info) {
	exonsPerTX = unlist(lapply(tx[subjectHits(ov)], length))
	exonsInd = data.frame(exonsPerTX, endInd = cumsum(exonsPerTX))
	exonsInd$startInd = exonsInd$endInd - exonsInd$exonsPerTX + 1
	for (trx in 1:nrow(exonsInd)) {
		for (row in exonsInd$startInd[trx]:exonsInd$endInd[trx]) {
			polygon(info$x[row, ], info$y[row, ], lwd = 1.2, col = gg_colors[trx])
		}
	}
}

pdf("isoform_structures_top_shifts.pdf",h=4,w=12)
par(mar=c(6,15,4,6))

for (i in 6:9) {
gene = top6[i]

## coordinates and TXs of gene
geneInfo = geneMap[gene,]
main = geneInfo$Symbol
genei = paste0(rownames(geneInfo),"\n",geneInfo$Symbol)
txsOfGene = unlist(geneInfo$gencodeTx)
chrom = geneInfo$seqnames
minRange = geneInfo$start
maxRange = geneInfo$end
strnd=geneInfo$strand
gr = GRanges(chrom, IRanges(minRange,maxRange), strnd)
yLabelCoords = -1:-length(txsOfGene)  ## negative y-coords so first tx is on top

## plot
    ov <- findOverlaps(gr, tx)
    txList <- split(tx[subjectHits(ov)], queryHits(ov))
	txIDs <- names(tx[subjectHits(ov)])
    poly.data <- lapply(txList, derfinderPlot:::.plotData)[[1]]
	if(strnd=="+") { 
		poly.data$intron$y = - poly.data$intron$y
		poly.data$exon$y = - poly.data$exon$y	}  ## negative so first is on top 
    yrange <- range(poly.data$exon$y)
    xrange <- c(start(gr),end(gr))
    plot(0,0, type="n", main=genei, xlim = xrange , ylim = yrange + c(-0.75, 0.75),
        yaxt="n", ylab="", xlab=as.character(seqnames(gr)),
        cex.axis = 1.2, cex.lab = 1.4, cex.main = 1.5)
	axis(2, at=yLabelCoords, labels=txIDs, las=1, cex.axis = 1.2)
    derfinderPlot:::.plotPolygon(poly.data$intron, 'gray90')
	.plotPolygonColorsExon(poly.data$exon)
}

dev.off()






##################################
######     original ggplot version
##################################
## barplot by genome and condition
library(ggplot2)
library(gridExtra)

gg_colors = c("#a38692","#70b670","#3290d8","#be9f52","#b12a2d","#7c52be","#df9c7f","#5e756b","#bbbed2",
		"#b8c85e","#7a633e","#fb895a","#4d5382")

## number to plot:
n=100

##
majorIsoAmountP2 = majorIsoAmountP2[order(rowMeans(majorIsoAmountP2),decreasing=TRUE),]
top6 = rownames(majorIsoAmountP2)[1:n]
top6tpm = txTpms2[top6]
top6percent = txPercents2[top6]
top6stats = gVoomStats[top6,]


pdf("isoform_top_shifts_top100.pdf", h=6, w=15, onefile=TRUE)
layout(matrix(c(1,1,1,1, 2,2,2), nc= 1))
for (i in 1:n) {
genei = paste0(rownames(top6stats)[i],"\n",top6stats$Symbol[i])
percenti = top6percent[[i]]
tpmi = top6tpm[[i]]
nTx = nrow(percenti)

means = data.frame(txID=rep(rownames(percenti),4),
				cond=c(rep("Accelerated\nDorsal",nTx),rep("NPC",nTx),rep("Rosette",nTx),rep("Neuron",nTx)),
				rate = c(percenti))
means$cond = factor(means$cond, levels=c("Accelerated\nDorsal","NPC","Rosette","Neuron") )
means[which(means==0, arr.ind=TRUE)] = 0.005 ## so it shows up on plot

meansTPM = data.frame(txID=rep(rownames(tpmi),4),
				cond=c(rep("Accelerated\nDorsal",nTx),rep("NPC",nTx),rep("Rosette",nTx),rep("Neuron",nTx)),
				rate = log2(unlist(c(tpmi+1))) )
meansTPM$cond = factor(meansTPM$cond, levels=c("Accelerated\nDorsal","NPC","Rosette","Neuron") )
meansTPM[which(meansTPM==0, arr.ind=TRUE)] = 0.005 ## so it shows up on plot


g1 = ggplot(means, aes(x=cond, y=rate, fill=factor(txID)))+
  geom_bar(stat="identity", position="dodge")+
  scale_fill_manual(values=gg_colors[1:nrow(percenti)], name="Transcript", 
		labels=rownames(percenti))+
  ggtitle(genei)+
  xlab("")+
  ylab("Relative proportion expressed")+
  ylim(0,1)+
  theme(panel.grid.major.x = element_blank(),
		plot.title=element_text(size=16, face="bold",color=1,hjust=.5),
		axis.text=element_text(size=14, face="bold",color=1),
		axis.title=element_text(size=18, face="bold",color=1),
		legend.text=element_text(size=12), 
		legend.title=element_text(size=13, face="bold"),
		#legend.position=c(.4,.92),
		legend.background = element_rect(fill=alpha('gray', 0.0)),
		legend.box.background = element_rect(color="darkgray"))

g2 = ggplot(meansTPM, aes(x=cond, y=rate, fill=factor(txID)))+
  geom_bar(stat="identity", position="dodge")+
  scale_fill_manual(values=gg_colors[1:nrow(tpmi)], name="Transcript", 
		labels=rownames(tpmi))+
  ggtitle(genei)+
  xlab("")+
  ylab("log2(TPM+1)")+
  #ylim(0,1)+
  theme(panel.grid.major.x = element_blank(),
		plot.title=element_text(size=16, face="bold",color=1,hjust=.5),
		axis.text=element_text(size=14, face="bold",color=1),
		axis.title=element_text(size=18, face="bold",color=1),
		legend.text=element_text(size=12), 
		legend.title=element_text(size=13, face="bold"),
		#legend.position=c(.4,.92),
		legend.background = element_rect(fill=alpha('gray', 0.0)),
		legend.box.background = element_rect(color="darkgray"))
g = grid.arrange(g2, g1, ncol=2)		
print(g)

}
dev.off()


