### load packages
library(jaffelab)
library(GenomicRanges)
library(rtracklayer)
library(SummarizedExperiment)
library(RColorBrewer)
library(rafalib)
library(edgeR)
library(limma)

MAINDIR = '/dcl01/lieber/ajaffe/lab/libd_stem_timecourse'
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n146.rda"))

geneCounts = assays(rse_gene)$counts
geneRpkm = getRPKM(rse_gene)
geneMap = rowData(rse_gene)

pd = colData(rse_gene)
pd$REP = as.numeric(pd$BIO_REP)
pd$REP[is.na(pd$REP)] = 1
pd$DAY = ordered(pd$DAY, levels=sort(as.numeric(levels(as.factor(pd$DAY)))) )
pd$DIV = as.numeric(as.character(pd$DAY))
pd$EXPERIMENT[pd$EXPERIMENT=="TIMECOURSE1"] = "TIMECOURSE"

## don't use technical replicates or RENEW controls in analyses
bioInd = which(pd$Class == "Naked genomes" & pd$CONDITION!="RENEW")
pd = pd[bioInd,]
geneRpkm = geneRpkm[,bioInd]
geneCounts = geneCounts[,bioInd]

pd$COND = pd$CONDITION
pd$COND[grep("NEURON",pd$COND)] = "NEURON"  ## combine on/off astros
pd$COND = as.factor(pd$COND)
pd$COND = factor(pd$COND,levels(pd$COND)[c(1,3,4,2)])  ## put levels in order

###########################
#### filter expression ####
geneMap$meanExprs = rowMeans(geneRpkm)
gIndex = which(geneMap$meanExprs > 0.1)
geneRpkm = geneRpkm[gIndex,]
geneCounts = geneCounts[gIndex,]
geneMap = geneMap[gIndex,]
map = geneMap


############################
## do voom analysis ########

##### fix model
## using ordered factor
mod = model.matrix(~COND + factor(LINE) + totalAssignedGene, data=pd)
mod0 = model.matrix(~factor(LINE) + totalAssignedGene, data=pd)

##############
## gene 
dge = DGEList(counts = geneCounts, genes = as.data.frame(geneMap))
dge = calcNormFactors(dge)
v = voom(dge, mod, plot=FALSE)
fit = lmFit(v)
ebGene = ebayes(fit)	
ffGene = topTable(eBayes(fit), coef=2:4, n = nrow(dge))
ffGene$bonf.P.Val = p.adjust(ffGene$P.Value, "bonf")
ffGene = ffGene[rownames(ebGene$t),]

voomStats = cbind(ffGene, ebGene$t[,2:4], ebGene$p.value[,2:4])
names(voomStats)[19:21] = paste0("t_", levels(pd$COND)[-1])
names(voomStats)[22:24] = paste0("pvalue_", levels(pd$COND)[-1])

save(voomStats, file="rda/voom_factor.rda")

###

############################
## plot top voomStats

# Colors by Donor
pal = c(colorRampPalette(brewer.pal(10,"Paired")[1:2])(4), 
		colorRampPalette(brewer.pal(10,"Paired")[3:4])(3), 
		brewer.pal(10,"Paired")[9:10], brewer.pal(10,"Paired")[5:8])
pd$lineCol = pal[as.numeric(factor(pd$LINE))]


sigOrderMat = as.data.frame(apply(voomStats[,22:24], 2, 
                                  function(x) order(x)[1:200]))
yExprs = as.matrix(log2(geneRpkm+1))

#### NPC ####
ooL = sigOrderMat$pvalue_NPC
pdf("ordered_effects_NPC.pdf",h=6,w=8)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
for(i in ooL) {
 plot(as.numeric(pd$DAY),yExprs[i,] , xaxt="n",
	   pch = 21, bg=pd$lineCol,
       cex=2,xlab="Day",
       ylab="log2(Exprs + 1)",
       main = paste0(geneMap$Symbol[i], "\n", geneMap$gencodeID[i]) )
 axis(1, at=1:1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
 legend("top", paste0("p=",signif(voomStats$pvalue_NPC[i],3)))
 if (i==ooL[1]) { legend("bottomleft", levels(factor(pd$LINE)),
	pch = 15, col = pal, cex=1) }
}
dev.off()

#### ROSETTE ####
ooL = sigOrderMat$pvalue_ROSETTE
pdf("ordered_effects_ROSETTE.pdf",h=6,w=8)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
for(i in ooL) {
 plot(as.numeric(pd$DAY),yExprs[i,] , xaxt="n",
	   pch = 21, bg=pd$lineCol,
       cex=2,xlab="Day",
       ylab="log2(Exprs + 1)",
       main = paste0(geneMap$Symbol[i], "\n", geneMap$gencodeID[i]) )
 axis(1, at=1:1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
 legend("top", paste0("p=",signif(voomStats$pvalue_ROSETTE[i],3)))
 if (i==ooL[1]) { legend("bottomleft", levels(factor(pd$LINE)),
	pch = 15, col = pal, cex=1) }
}
dev.off()

#### NEURON ####
ooL = sigOrderMat$pvalue_NEURON
pdf("ordered_effects_NEURON.pdf",h=6,w=8)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
for(i in ooL) {
 plot(as.numeric(pd$DAY),yExprs[i,] , xaxt="n",
	   pch = 21, bg=pd$lineCol,
       cex=2,xlab="Day",
       ylab="log2(Exprs + 1)",
       main = paste0(geneMap$Symbol[i], "\n", geneMap$gencodeID[i]) )
 axis(1, at=1:1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
 legend("top", paste0("p=",signif(voomStats$pvalue_NEURON[i],3)))
 if (i==ooL[1]) { legend("bottomleft", levels(factor(pd$LINE)),
	pch = 15, col = pal, cex=1) }
}
dev.off()







