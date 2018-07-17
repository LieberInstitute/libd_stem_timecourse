### load packages
library(jaffelab)
library(GenomicRanges)
library(rtracklayer)
library(SummarizedExperiment)
library(RColorBrewer)
library(rafalib)
library(edgeR)
library(limma)

#############################
#### load feature counts ####
MAINDIR = '/dcl01/lieber/ajaffe/lab/libd_stem_timecourse'
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n146.rda"))
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseExon_n146.rda"))
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseJxn_n146.rda"))
load(file.path(MAINDIR,"data/libd_stemcell_rseTx_counts.rda"))

# #########################
# #### add transcripts #### ## need counts for voom
# txPath = "/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/Salmon_tx/"
# txFiles = paste0(txPath, rse_gene$SAMPLE_ID, "/quant.sf")
# names(txFiles) = rse_gene$SampleID
# txList = lapply(txFiles, read.delim, row.names=1, as.is=TRUE)
# txTpm = sapply(txList, "[[", "TPM")
# txCounts = sapply(txList, "[[", "NumReads")
# ##get names of transcripts
# txNames = rownames(txList[[1]])
# txNames = as.character(txNames)
# txMap = t(ss(txNames, "\\|",c(1,7,2,6,8)))
# txMap = as.data.frame(txMap)
# colnames(txMap) = c("gencodeTx","Length","gencodeID","Symbol","gene_type")
# rownames(txMap) = rownames(txTpm) = rownames(txCounts) = txMap$gencodeTx
# ## make tx RSE
# rse_tx = SummarizedExperiment(
	# assays = list('tpm'=txTpm, 'counts'=txCounts),
    # colData = as.data.frame(colData(rse_gene)), 
	# rowData = txMap)
# save(rse_tx, file="/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/data/libd_stemcell_rseTx_counts.rda")

###################
#### functions ####
## recount getRPKM version
getRPKM <- function(rse, length_var = 'bp_length', mapped_var = NULL) {
    mapped <- if(!is.null(mapped_var)) colData(rse)[, mapped_var] else colSums(assays(rse)$counts)      
    bg <- matrix(mapped, ncol = ncol(rse), nrow = nrow(rse), byrow = TRUE)   
    len <- if(!is.null(length_var)) rowData(rse)[, length_var] else width(rowRanges(rse))   
    wid <- matrix(len, nrow = nrow(rse), ncol = ncol(rse), byrow = FALSE)
    assays(rse)$counts / (wid/1000) / (bg/1e6)
}
getRPM = function(rse, target = 80e6) {
	require(SummarizedExperiment)
	mapped <- colSums(assays(rse)$counts) 
	bg = matrix(rep(mapped/target), nc = ncol(rse), nr = nrow(rse),	byrow=TRUE)
	assays(rse)$counts/bg
}

########################
#### subset samples ####
dropInd = which(rse_gene$Class != "Naked genomes" | rse_gene$CONDITION %in% c("RENEW","NEURONS_ALONE"))

rse_gene = rse_gene[,-dropInd]  ## n=106
rse_exon = rse_exon[,-dropInd]
rse_jxn = rse_jxn[,-dropInd]
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
geneRpkm = getRPKM(rse_gene, length_var="Length")
exonRpkm = getRPKM(rse_exon, length_var="Length")
jRpkm = getRPM(rse_jxn, target = 10e6)
txTpm = assays(rse_tx)$tpm

mcols(rse_gene)$meanExprs = rowMeans(geneRpkm)
rse_gene = rse_gene[which(mcols(rse_gene)$meanExprs > 0.1),]

mcols(rse_exon)$meanExprs = rowMeans(exonRpkm)
rse_exon = rse_exon[which(mcols(rse_exon)$meanExprs > 0.3),]

mcols(rse_jxn)$meanExprs = rowMeans(jRpkm)
rse_jxn = rse_jxn[which(mcols(rse_jxn)$meanExprs > 0.3),]

mcols(rse_tx)$meanExprs = rowMeans(txTpm)
rse_tx = rse_tx[which(mcols(rse_tx)$meanExprs > 0.3),]


########################
#### counts objects ####
geneCounts = assays(rse_gene)$counts
geneMap = as.data.frame(rowRanges(rse_gene))
exonCounts = assays(rse_exon)$counts
exonMap = as.data.frame(rowRanges(rse_exon))
jCounts = assays(rse_jxn)$counts
jMap = as.data.frame(rowRanges(rse_jxn))
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



# ############################################
# ############ voom - genes ##################
# ############################################

# # using contrasts
# mod = model.matrix(~0 + COND + LINE + totalAssignedGene, data=pd)
# colnames( mod ) <- c(levels(pd$COND), colnames(mod)[-seq(length(levels(pd$COND)))])
# cmtx <- makeContrasts( "NPC-ACC_DORSAL", "ROSETTE-NPC",'NEURON-ROSETTE', levels= mod)

# dgeGene = DGEList(counts = geneCounts, genes = as.data.frame(geneMap))
# dgeGene = calcNormFactors(dgeGene)
# vGene = voom(dgeGene, mod, plot=FALSE)
# fitGene = lmFit(vGene)
# ebGene = ebayes(contrasts.fit(fitGene,cmtx))
# ffGene = topTable(eBayes(contrasts.fit(fitGene,cmtx)), coef=1:3, n = nrow(dgeGene)) #by condition
# ffGene$bonf.P.Val = p.adjust(ffGene$P.Value, "bonf")
# ffGene = ffGene[rownames(ebGene$t),]

# significance levels
# pvalMat = as.matrix(ebGene$p.value)
# qvalMat = pvalMat
# qvalMat[,1:3] = p.adjust(pvalMat[,1:3],method="fdr") 
# colnames(pvalMat) = paste0("p_",colnames(pvalMat))
# colnames(qvalMat) = paste0("q_",colnames(qvalMat))

# sig = cbind(pvalMat,qvalMat)

# total genes: 
# nrow(geneCounts)
# [1] 25,466

# table(ffGene$"adj.P.Val" < 0.01)
# FALSE  TRUE
 # 5246 20220

# colSums(sig<0.01)
# q_NPC-ACC_DORSAL   q_ROSETTE-NPC  q_NEURON-ROSETTE
# 9067              1994             12951

# gVoomStats = cbind(ffGene, ebGene$t, pvalMat, qvalMat)
# names(gVoomStats)[24:26] = paste0("t_", names(gVoomStats)[24:26])


### Percentage of variance
library(variancePartition)
library(Matrix)

form <- ~ COND + Donor + LINE + totalAssignedGene
elist = vGene$E

varPart <- fitVarPartModel(elist, form, as.data.frame(pd), showWarnings = FALSE)

varGene = extractVarPart(varPart, showWarnings = FALSE)
apply(varGene, 2, mean)

apply(varGene, 2, sd)

p = plotVarPart(varGene, label.angle = 60)

pdf("varianceExplained_gene2.pdf", h=6, w=8)
p + theme_bw(base_size=18) + theme(legend.position="none")
dev.off()

varGene$donorPLUSline = varGene$Donor + varGene$LINE
table(varGene$donorPLUSline > varGene$COND)
mean(varGene$donorPLUSline > varGene$COND)   # 0.2352548
FALSE  TRUE
19475  5991
table(varGene$Donor > varGene$COND)
mean(varGene$Donor > varGene$COND)			# 0.09915181
FALSE  TRUE
22941  2525
table(varGene$LINE > varGene$COND)
mean(varGene$LINE > varGene$COND)			# 0.1543234
FALSE  TRUE
21536  3930


# ##############################################
# ############## voom - exons ##################
# ##############################################

# ## same mod

# dge = DGEList(counts = exonCounts, genes = as.data.frame(exonMap))
# dge = calcNormFactors(dge)
# v = voom(dge, mod, plot=FALSE)
# fit = lmFit(v)
# eb = ebayes(contrasts.fit(fit,cmtx))
# ff = topTable(eBayes(contrasts.fit(fit,cmtx)), coef=1:3, n = nrow(dge)) #by condition
# ff$bonf.P.Val = p.adjust(ff$P.Value, "bonf")
# ff = ff[rownames(eb$t),]

# ## significance levels
# pvalMat = as.matrix(eb$p.value)
# qvalMat = pvalMat
# qvalMat[,1:3] = p.adjust(pvalMat[,1:3],method="fdr") 
# colnames(pvalMat) = paste0("p_",colnames(pvalMat))
# colnames(qvalMat) = paste0("q_",colnames(qvalMat))

# sig = cbind(pvalMat,qvalMat)

# ## total exons: 
# nrow(exonCounts)
# # [1] 356,663

# table(ff$"adj.P.Val" < 0.01)
# # FALSE   TRUE
# # 62,537 294,126

# colSums(sig<0.01)
# # q_NPC-ACC_DORSAL     q_ROSETTE-NPC 		q_NEURON-ROSETTE
# # 115,449              20,569           	189,437

# eVoomStats = cbind(ff, eb$t, pvalMat, qvalMat)
# names(eVoomStats)[24:26] = paste0("t_", names(eVoomStats)[24:26])



### Percentage of variance
library(variancePartition)

form <- ~ COND + LINE + totalAssignedGene
elist = v$E

varPart <- fitVarPartModel(elist, form, as.data.frame(pd), showWarnings = FALSE)

varExon = extractVarPart(varPart, showWarnings = FALSE)
apply(varExon, 2, mean)
             # COND              LINE totalAssignedGene         Residuals
       # 0.43095394        0.12443977        0.04358094        0.40102534
apply(varExon, 2, sd)
             # COND              LINE totalAssignedGene         Residuals
       # 0.28018158        0.10297282        0.05392823        0.21811132

pdf("varianceExplained_exon.pdf")
plotVarPart(varExon)
dev.off()



# ##############################################
# ############## voom - jxns ##################
# ##############################################

# ## same mod

# dge = DGEList(counts = jCounts, genes = as.data.frame(jMap))
# dge = calcNormFactors(dge)
# v = voom(dge, mod, plot=FALSE)
# fit = lmFit(v)
# eb = ebayes(contrasts.fit(fit,cmtx))
# ff = topTable(eBayes(contrasts.fit(fit,cmtx)), coef=1:3, n = nrow(dge)) #by condition
# ff$bonf.P.Val = p.adjust(ff$P.Value, "bonf")
# ff = ff[rownames(eb$t),]

# ## significance levels
# pvalMat = as.matrix(eb$p.value)
# qvalMat = pvalMat
# qvalMat[,1:3] = p.adjust(pvalMat[,1:3],method="fdr") 
# colnames(pvalMat) = paste0("p_",colnames(pvalMat))
# colnames(qvalMat) = paste0("q_",colnames(qvalMat))

# sig = cbind(pvalMat,qvalMat)

# ## total jxns: 
# nrow(jCounts)
# # [1] 281,497

# table(ff$"adj.P.Val" < 0.01)
# # FALSE   TRUE
 # # 99,840  181,657

# colSums(sig<0.01)
# # q_NPC-ACC_DORSAL   q_ROSETTE-NPC  	q_NEURON-ROSETTE
# # 57,074             9,816             107,110

# jVoomStats = cbind(ff, eb$t, pvalMat, qvalMat)
# names(jVoomStats)[30:32] = paste0("t_", names(jVoomStats)[30:32])



# ##############################################
# ############## voom - tx ##################
# ##############################################

# ## same mod

# dge = DGEList(counts = txCounts, genes = as.data.frame(txMap))
# dge = calcNormFactors(dge)
# v = voom(dge, mod, plot=FALSE)
# fit = lmFit(v)
# eb = ebayes(contrasts.fit(fit,cmtx))
# ff = topTable(eBayes(contrasts.fit(fit,cmtx)), coef=1:3, n = nrow(dge)) #by condition
# ff$bonf.P.Val = p.adjust(ff$P.Value, "bonf")
# ff = ff[rownames(eb$t),]

# ## significance levels
# pvalMat = as.matrix(eb$p.value)
# qvalMat = pvalMat
# qvalMat[,1:3] = p.adjust(pvalMat[,1:3],method="fdr") 
# colnames(pvalMat) = paste0("p_",colnames(pvalMat))
# colnames(qvalMat) = paste0("q_",colnames(qvalMat))

# sig = cbind(pvalMat,qvalMat)

# ## total transcripts: 
# nrow(txCounts)
# # [1] 73,140

# table(ff$"adj.P.Val" < 0.01)
# # FALSE   TRUE
 # # 22,407  50,733

# colSums(sig<0.01)
# # q_NPC-ACC_DORSAL    q_ROSETTE-NPC  	q_NEURON-ROSETTE
# # 12,616              2109             31,672

# txVoomStats = cbind(ff, eb$t, pvalMat, qvalMat)
# names(txVoomStats)[15:17] = paste0("t_", names(txVoomStats)[15:17])


# save(gVoomStats, eVoomStats, jVoomStats, txVoomStats, file="voomStats_4features.rda")



load("voomStats_4features.rda", verbose=TRUE)


###############
#### plots ####
###############

# Colors by Donor
pal = c(colorRampPalette(brewer.pal(10,"Paired")[1:2])(4), 
		colorRampPalette(brewer.pal(10,"Paired")[3:4])(3), 
		brewer.pal(10,"Paired")[9:10], brewer.pal(10,"Paired")[5:8])
pd$lineCol = pal[as.numeric(factor(pd$LINE))]

geneRpkm = getRPKM(rse_gene, length_var="Length")
exonRpkm = getRPKM(rse_exon, length_var="Length")
jRpkm = getRPM(rse_jxn, target = 10e6)
txTpm = assays(rse_tx)$tpm
gExprs = as.matrix(log2(geneRpkm+1))
eExprs = as.matrix(log2(exonRpkm+1))
jExprs = as.matrix(log2(jRpkm+1))
txExprs = as.matrix(log2(txTpm+1))



##############
#### gene ####
sigOrderMat = as.data.frame(apply(gVoomStats[,c(27:29,21)], 2, 
                                  function(x) order(x)[1:200]))

gPlotTop = function(ind=4) {
fileLab = c("1_AD_NPC","2_NPC_ROSETTE","3_ROSETTE_NEURON","4_overall")[ind]
voomStatsInd = c(27:29,21)[ind]
ooL = sigOrderMat[,ind]
pdf(paste0("gene_ordered_effects_contrasts_",fileLab,".pdf"),h=6,w=8)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
for(i in ooL) {
 plot(as.numeric(pd$DAY),gExprs[i,] , xaxt="n", pch = 21, 
	   bg=pd$lineCol, cex=2,xlab="Day", ylab="log2(Exprs + 1)",
       main = paste0(geneMap$Symbol[i], "\n", geneMap$gencodeID[i]) )
 axis(1, at=1:1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
 abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
 legend("top", paste0("p=",signif(gVoomStats[,voomStatsInd][i],3)), bg="white")
 if (i==ooL[1]) { legend("topright", levels(factor(pd$LINE)),
	pch = 15, col = pal, cex=1) }
}
dev.off()
}

for (j in 1:4) { gPlotTop(ind=j) }



##############
#### exon ####
sigOrderMat = as.data.frame(apply(eVoomStats[,c(27:29,21)], 2, 
                                  function(x) order(x)[1:200]))

ePlotTop = function(ind=4) {
fileLab = c("1_AD_NPC","2_NPC_ROSETTE","3_ROSETTE_NEURON","4_overall")[ind]
voomStatsInd = c(27:29,21)[ind]
ooL = sigOrderMat[,ind]
pdf(paste0("exon_ordered_effects_contrasts_",fileLab,".pdf"),h=6,w=8)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
for(i in ooL) {
 plot(as.numeric(pd$DAY),eExprs[i,] , xaxt="n", pch = 21, 
	   bg=pd$lineCol, cex=2,xlab="Day", ylab="log2(Exprs + 1)",
       main = paste0(exonMap$Symbol[i], " (exon) \n", 
	   exonMap$seqnames[i],":",exonMap$start[i],"-",exonMap$end[i],"(",exonMap$strand[i],")") )
 axis(1, at=1:1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
 abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
 legend("top", paste0("p=",signif(eVoomStats[,voomStatsInd][i],3)), bg="white")
 if (i==ooL[1]) { legend("topright", levels(factor(pd$LINE)),
	pch = 15, col = pal, cex=1) }
}
dev.off()
}

for (j in 1:4) { ePlotTop(ind=j) }



##############
#### jxn #####
sigOrderMat = as.data.frame(apply(jVoomStats[,c(33:35,27)], 2, 
                                  function(x) order(x)[1:200]))
jMap$Symbol[jMap$Class!="InGen"] = paste0("[",jMap$Class[jMap$Class!="InGen"],"]")
								  
jPlotTop = function(ind=4) {
fileLab = c("1_AD_NPC","2_NPC_ROSETTE","3_ROSETTE_NEURON","4_overall")[ind]
voomStatsInd = c(33:35,27)[ind]
ooL = sigOrderMat[,ind]
pdf(paste0("jxn_ordered_effects_contrasts_",fileLab,".pdf"),h=6,w=8)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
for(i in ooL) {
 plot(as.numeric(pd$DAY),jExprs[i,] , xaxt="n", pch = 21, 
	   bg=pd$lineCol, cex=2,xlab="Day", ylab="log2(Exprs + 1)",
       main = paste0(jMap$Symbol[i], " (junction) \n", 
	   jMap$seqnames[i],":",jMap$start[i],"-",jMap$end[i],"(",jMap$strand[i],")") )
 axis(1, at=1:1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
 abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
 legend("top", paste0("p=",signif(jVoomStats[,voomStatsInd][i],3)), bg="white")
 if (i==ooL[1]) { legend("topright", levels(factor(pd$LINE)),
	pch = 15, col = pal, cex=1) }
}
dev.off()
}

for (j in 1:4) { jPlotTop(ind=j) }




##############
#### tx #####
sigOrderMat = as.data.frame(apply(txVoomStats[,c(18:20,12)], 2, 
                                  function(x) order(x)[1:200]))

								  
txPlotTop = function(ind=4) {
fileLab = c("1_AD_NPC","2_NPC_ROSETTE","3_ROSETTE_NEURON","4_overall")[ind]
voomStatsInd = c(18:20,12)[ind]
ooL = sigOrderMat[,ind]
pdf(paste0("tx_ordered_effects_contrasts_",fileLab,".pdf"),h=6,w=8)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
for(i in ooL) {
 plot(as.numeric(pd$DAY),txExprs[i,] , xaxt="n", pch = 21, 
	   bg=pd$lineCol, cex=2,xlab="Day", ylab="log2(Exprs + 1)",
       main = paste0(txMap$Symbol[i], "\n", txMap$gencodeTx[i]) )
 axis(1, at=1:1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
 abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
 legend("top", paste0("p=",signif(txVoomStats[,voomStatsInd][i],3)), bg="white")
 if (i==ooL[1]) { legend("topright", levels(factor(pd$LINE)),
	pch = 15, col = pal, cex=1) }
}
dev.off()
}

for (j in 1:4) { txPlotTop(ind=j) }






############################
## plot top voomStatsC

load("voomStats_4features.rda", verbose=TRUE)


gVoomStats$AD_meanExprs = rowMeans(geneRpkm[,which(pd$CONDITION=="ACC_DORSAL(2)")] )
gVoomStats$NPC_meanExprs = rowMeans(geneRpkm[,which(pd$CONDITION=="NPC")] )
gVoomStats$ROSE_meanExprs = rowMeans(geneRpkm[,which(pd$CONDITION=="ROSETTE")] )
gVoomStats$NEUR_meanExprs = rowMeans(geneRpkm[,which(pd$CONDITION=="NEURONS_PLUS_ASTROS")] )
gVoomStats = gVoomStats[, -c(7:8,11:12,15:19)]
#### AD-NPC
gSig1 = gVoomStats[order(gVoomStats$"p_NPC-ACC_DORSAL"),]
gSig1 = gSig1[gSig1$"q_NPC-ACC_DORSAL" < 0.01, ]
write.csv(gSig1, file="csv/gene_voom_FDR01_1_AD_NPC.csv", quote=FALSE)
#### NPC-ROSETTE
gSig2 = gVoomStats[order(gVoomStats$"p_ROSETTE-NPC"),]
gSig2 = gSig2[gSig2$"q_ROSETTE-NPC" < 0.01, ]
write.csv(gSig2, file="csv/gene_voom_FDR01_2_NPC_ROSETTE.csv", quote=FALSE)
#### ROSETTE-NEURONS
gSig3 = gVoomStats[order(gVoomStats$"p_NEURON-ROSETTE"),]
gSig3 = gSig3[gSig3$"q_NEURON-ROSETTE" < 0.01, ]
write.csv(gSig3, file="csv/gene_voom_FDR01_3_ROSETTE_NEURON.csv", quote=FALSE)
#### OVERALL
gSig4 = gVoomStats[order(gVoomStats$"P.Value"),]
gSig4 = gSig4[gSig4$"adj.P.Val" < 0.01, ]
write.csv(gSig4, file="csv/gene_voom_FDR01_4_OVERALL.csv", quote=FALSE)



eVoomStats$AD_meanExprs = rowMeans(exonRpkm[,which(pd$CONDITION=="ACC_DORSAL(2)")] )
eVoomStats$NPC_meanExprs = rowMeans(exonRpkm[,which(pd$CONDITION=="NPC")] )
eVoomStats$ROSE_meanExprs = rowMeans(exonRpkm[,which(pd$CONDITION=="ROSETTE")] )
eVoomStats$NEUR_meanExprs = rowMeans(exonRpkm[,which(pd$CONDITION=="NEURONS_PLUS_ASTROS")] )
eVoomStats = eVoomStats[, -c(8,11:12,15:19)]
#### AD-NPC
eSig1 = eVoomStats[order(eVoomStats$"p_NPC-ACC_DORSAL"),]
eSig1 = eSig1[eSig1$"q_NPC-ACC_DORSAL" < 0.01, ]
write.csv(eSig1, file="csv/exon_voom_FDR01_1_AD_NPC.csv", quote=FALSE)
#### NPC-ROSETTE
eSig2 = eVoomStats[order(eVoomStats$"p_ROSETTE-NPC"),]
eSig2 = eSig2[eSig2$"q_ROSETTE-NPC" < 0.01, ]
write.csv(eSig2, file="csv/exon_voom_FDR01_2_NPC_ROSETTE.csv", quote=FALSE)
#### ROSETTE-NEURONS
eSig3 = eVoomStats[order(eVoomStats$"p_NEURON-ROSETTE"),]
eSig3 = eSig3[eSig3$"q_NEURON-ROSETTE" < 0.01, ]
write.csv(eSig3, file="csv/exon_voom_FDR01_3_ROSETTE_NEURON.csv", quote=FALSE)
#### OVERALL
eSig4 = eVoomStats[order(eVoomStats$"P.Value"),]
eSig4 = eSig4[eSig4$"adj.P.Val" < 0.01, ]
write.csv(eSig4, file="csv/exon_voom_FDR01_4_OVERALL.csv", quote=FALSE)


jVoomStats$AD_meanExprs = rowMeans(jRpkm[,which(pd$CONDITION=="ACC_DORSAL(2)")] )
jVoomStats$NPC_meanExprs = rowMeans(jRpkm[,which(pd$CONDITION=="NPC")] )
jVoomStats$ROSE_meanExprs = rowMeans(jRpkm[,which(pd$CONDITION=="ROSETTE")] )
jVoomStats$NEUR_meanExprs = rowMeans(jRpkm[,which(pd$CONDITION=="NEURONS_PLUS_ASTROS")] )
jVoomStats = jVoomStats[, -c(9:13,22:25)]
#### AD-NPC
jSig1 = jVoomStats[order(jVoomStats$"p_NPC-ACC_DORSAL"),]
jSig1 = jSig1[jSig1$"q_NPC-ACC_DORSAL" < 0.01, ]
write.csv(jSig1, file="csv/jxn_voom_FDR01_1_AD_NPC.csv", quote=FALSE)
#### NPC-ROSETTE
jSig2 = jVoomStats[order(jVoomStats$"p_ROSETTE-NPC"),]
jSig2 = jSig2[jSig2$"q_ROSETTE-NPC" < 0.01, ]
write.csv(jSig2, file="csv/jxn_voom_FDR01_2_NPC_ROSETTE.csv", quote=FALSE)
#### ROSETTE-NEURONS
jSig3 = jVoomStats[order(jVoomStats$"p_NEURON-ROSETTE"),]
jSig3 = jSig3[jSig3$"q_NEURON-ROSETTE" < 0.01, ]
write.csv(jSig3, file="csv/jxn_voom_FDR01_3_ROSETTE_NEURON.csv", quote=FALSE)
#### OVERALL
jSig4 = jVoomStats[order(jVoomStats$"P.Value"),]
jSig4 = jSig4[jSig4$"adj.P.Val" <= 0.01, ]
write.csv(jSig4, file="csv/jxn_voom_FDR01_4_OVERALL.csv", quote=FALSE)


txVoomStats$AD_meanExprs = rowMeans(txTpm[,which(pd$CONDITION=="ACC_DORSAL(2)")] )
txVoomStats$NPC_meanExprs = rowMeans(txTpm[,which(pd$CONDITION=="NPC")] )
txVoomStats$ROSE_meanExprs = rowMeans(txTpm[,which(pd$CONDITION=="ROSETTE")] )
txVoomStats$NEUR_meanExprs = rowMeans(txTpm[,which(pd$CONDITION=="NEURONS_PLUS_ASTROS")] )
txVoomStats = txVoomStats[, -c(1,7:10)]
#### AD-NPC
txSig1 = txVoomStats[order(txVoomStats$"p_NPC-ACC_DORSAL"),]
txSig1 = txSig1[txSig1$"q_NPC-ACC_DORSAL" < 0.01, ]
write.csv(txSig1, file="csv/tx_voom_FDR01_1_AD_NPC.csv", quote=FALSE)
#### NPC-ROSETTE
txSig2 = txVoomStats[order(txVoomStats$"p_ROSETTE-NPC"),]
txSig2 = txSig2[txSig2$"q_ROSETTE-NPC" < 0.01, ]
write.csv(txSig2, file="csv/tx_voom_FDR01_2_NPC_ROSETTE.csv", quote=FALSE)
#### ROSETTE-NEURONS
txSig3 = txVoomStats[order(txVoomStats$"p_NEURON-ROSETTE"),]
txSig3 = txSig3[txSig3$"q_NEURON-ROSETTE" < 0.01, ]
write.csv(txSig3, file="csv/tx_voom_FDR01_3_ROSETTE_NEURON.csv", quote=FALSE)
#### OVERALL
txSig4 = txVoomStats[order(txVoomStats$"P.Value"),]
txSig4 = txSig4[txSig4$"adj.P.Val" < 0.01, ]
write.csv(txSig4, file="csv/tx_voom_FDR01_4_OVERALL.csv", quote=FALSE)



### Table of numbers
sigFeatFDR01 = data.frame(AD_NPC=NA, NPC_ROSETTE=NA, ROSETTE_NEURON=NA, TIMECOURSE=NA)
sigFeatFDR01[1,] = c(nrow(gSig1),nrow(gSig2),nrow(gSig3),nrow(gSig4))
sigFeatFDR01[2,] = c(nrow(eSig1),nrow(eSig2),nrow(eSig3),nrow(eSig4))
sigFeatFDR01[3,] = c(nrow(jSig1),nrow(jSig2),nrow(jSig3),nrow(jSig4))
sigFeatFDR01[4,] = c(nrow(txSig1),nrow(txSig2),nrow(txSig3),nrow(txSig4))
rownames(sigFeatFDR01)=c("Genes","Exons","Junctions","Transcripts")
sigFeatFDR01$total_tested = c(nrow(geneCounts),nrow(exonCounts),nrow(jCounts),nrow(txCounts))

sigFeatFDR01
            # AD_NPC NPC_ROSETTE ROSETTE_NEURON TIMECOURSE total_tested
# Genes         9067        1994          12951      20220        25466
# Exons       115449       20569         189437     294126       356663
# Junctions    57074        9816         107110     181657       281497
# Transcripts  12616        2109          31672      50733        73140

sigGeneFDR01 = data.frame(AD_NPC=NA, NPC_ROSETTE=NA, ROSETTE_NEURON=NA, TIMECOURSE=NA)
sigGeneFDR01[1,] = c(nrow(gSig1),nrow(gSig2),nrow(gSig3),nrow(gSig4))
sigGeneFDR01[2,] = c(length(unique(eSig1$gencodeID)),length(unique(eSig2$gencodeID)),length(unique(eSig3$gencodeID)),length(unique(eSig4$gencodeID)))
sigGeneFDR01[3,] = c(length(unique(jSig1$newGeneID)),length(unique(jSig2$newGeneID)),length(unique(jSig3$newGeneID)),length(unique(jSig4$newGeneID)))
sigGeneFDR01[4,] = c(length(unique(txSig1$gencodeID)),length(unique(txSig2$gencodeID)),length(unique(txSig3$gencodeID)),length(unique(txSig4$gencodeID)))
rownames(sigGeneFDR01)=c("Genes","Exons","Junctions","Transcripts")

sigGeneFDR01
            # AD_NPC NPC_ROSETTE ROSETTE_NEURON TIMECOURSE
# Genes         9067        1994          12951      20220
# Exons        10221        3156          14173      17954
# Junctions     9665        3065          13185      16861
# Transcripts   7868        1752          14262      18359


sigJxnFDR01 = cbind(table(jSig1$Class),table(jSig2$Class),table(jSig3$Class),table(jSig4$Class),table(jMap$Class))
colnames(sigJxnFDR01) = c("AD_NPC","NPC_ROSETTE","ROSETTE_NEURON","TIMECOURSE","total_tested")
rownames(sigJxnFDR01) = names(table(jSig1$Class))
sigJxnFDR01 = sigJxnFDR01[c(3,1,2,4),]

sigJxnFDR01
            # AD_NPC NPC_ROSETTE ROSETTE_NEURON TIMECOURSE total_tested
# InGen        50106        8305          92625     151194       202339
# AltStartEnd   3325         704           7199      15002        38147
# ExonSkip      1798         334           3333       7298        18303
# Novel         1845         473           3953       8163        22708


write.csv(sigFeatFDR01, file="csv/summary_4features_FDR01.csv")
write.csv(sigGeneFDR01, file="csv/summary_4features_unique_genes_FDR01.csv")
write.csv(sigJxnFDR01, file="csv/summary_jxn_classes_FDR01.csv")





