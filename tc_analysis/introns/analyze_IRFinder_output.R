### load packages
library(jaffelab)
library(GenomicRanges)
library(rtracklayer)
library(SummarizedExperiment)
library(RColorBrewer)
library(rafalib)
library(edgeR)
library(limma)
library(BiocParallel)


######################
## load data
MAINDIR = '/dcl01/lieber/ajaffe/lab/libd_stem_timecourse'
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n146.rda"))

geneCounts = assays(rse_gene)$counts
geneMap = as.data.frame(rowRanges(rse_gene))

## recount getRPKM version
getRPKM <- function(rse, length_var = 'bp_length', mapped_var = NULL) {
    mapped <- if(!is.null(mapped_var)) colData(rse)[, mapped_var] else colSums(assays(rse)$counts)      
    bg <- matrix(mapped, ncol = ncol(rse), nrow = nrow(rse), byrow = TRUE)   
    len <- if(!is.null(length_var)) rowData(rse)[, length_var] else width(rowRanges(rse))   
    wid <- matrix(len, nrow = nrow(rse), ncol = ncol(rse), byrow = FALSE)
    assays(rse)$counts / (wid/1000) / (bg/1e6)
}
geneRpkm = getRPKM(rse_gene, length_var="Length")


pd = colData(rse_gene)
pd$DAY = ordered(pd$DAY, levels=sort(as.numeric(levels(as.factor(pd$DAY)))) )
pd$DIV = as.numeric(as.character(pd$DAY))
pd$EXPERIMENT[pd$EXPERIMENT=="TIMECOURSE1"] = "TIMECOURSE"
pd$LINE = gsub('-','_',pd$LINE)

## don't use technical replicates
bioInd = which(pd$Class == "Naked genomes")
pd = pd[bioInd,]
geneRpkm = geneRpkm[,bioInd]
geneCounts = geneCounts[,bioInd]   # n=128
## don't use RENEW controls or NEURONS_ALONE in analyses
dropInd = which(pd$CONDITION %in% c("RENEW","NEURONS_ALONE"))
pd = pd[-dropInd,]
geneRpkm = geneRpkm[,-dropInd]
geneCounts = geneCounts[,-dropInd]   # n=106

## clean CONDITION
pd$COND = as.factor(pd$CONDITION)
pd$COND = factor(pd$COND,levels(pd$COND)[c(1,3,4,2)])  ## put levels in order

# Colors by Donor
pal = c(colorRampPalette(brewer.pal(10,"Paired")[1:2])(4), 
		colorRampPalette(brewer.pal(10,"Paired")[3:4])(3), 
		brewer.pal(10,"Paired")[9:10], brewer.pal(10,"Paired")[5:8])
pd$lineCol = pal[as.numeric(factor(pd$LINE))]


######################
######################

intronDir = "/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/Counts/intron"
intronFn = file.path(intronDir, paste0(pd$SAMPLE_ID, '_Gencode.v25.hg38_Introns.counts.summary'))

# number of reads assigned to introns
intronStatList = lapply(intronFn, read.delim,row.names=1)
intronStats = do.call("cbind", intronStatList)
colnames(intronStats) = pd$SAMPLE_ID

pd$assignedIntron = as.numeric(intronStats[1,])
pd$assignedIntronRate = as.numeric(intronStats[1,] / colSums(intronStats))
pd$assignedIntronRate2 = as.numeric(intronStats[1,] / (intronStats[1,]+intronStats[2,]+intronStats[4,]))

## p-value of difference between npc/rosette and neuron intron rates
pd2 = pd[which(pd$CONDITION != "ACC_DORSAL(2)"),]
pd2$COND = ifelse(pd2$CONDITION=="NEURONS_PLUS_ASTROS", "NEURON", "NPCROSETTE")
fit = lm(assignedIntronRate2 ~ COND + LINE, data=pd2)
summary(fit)
summary(fit)$coef[2,4]

pdf("intron_gene.pdf", h=5,w=6)
plot(pd$totalAssignedGene, pd$assignedIntronRate2, col = pd$lineCol)

  boxplot(pd$totalAssignedGene ~ pd$DAY, xlab="Day", ylim=c(.1,.8),
		main="gene assignment rate", outline=FALSE)
  points(pd$totalAssignedGene ~ jitter(as.numeric(factor(pd$DAY)), 
		  amount=0.15), cex=1.5,
		  pch=21, bg = pd$lineCol)
		  
  boxplot(pd$assignedIntronRate2 ~ pd$DAY, xlab="Day", ylim=c(0.03,.22),
		main="Overall Intron Assigment Rate", outline=FALSE)
  points(pd$assignedIntronRate2 ~ jitter(as.numeric(factor(pd$DAY)), 
		  amount=0.15), cex=1.5,
		  pch=21, bg = pd$lineCol)		  

dev.off()

pdf("intron_assignment_over_time.pdf", h=5,w=6)
par(mar=c(5,4,3,2),cex.axis=1.3,cex.lab=1.5,cex.main=1.5)
  boxplot(pd$assignedIntronRate2 ~ pd$DAY, xlab="Day", ylim=c(0.03,.22),
		main="Overall Intron Assigment Rate", outline=FALSE)
  points(pd$assignedIntronRate2 ~ jitter(as.numeric(factor(pd$DAY)), 
		  amount=0.15), cex=1.5,
		  pch=21, bg = pd$lineCol)
  abline(v=c(4.5,5.5,6.5),col="gray",lty=2)
dev.off()

pdf("totalAssignedIntron.pdf", h=12,w=8)
par(mfrow=c(2,1))
boxplot(pd$assignedIntronRate ~ pd$DAY, main="Introns: Assigned/colSums(intronStats) \n including Unassigned_Multimapping")
abline(v=c(4.5,5.5,6.5),col="gray",lty=2)
boxplot(pd$assignedIntronRate2 ~ pd$DAY, main="Introns: Assigned/colSums(intronStats) \n not including Unassigned_Multimapping")
abline(v=c(4.5,5.5,6.5),col="gray",lty=2)
dev.off()

######################
######################


######################
## create tables

idir = "/dcl01/lieber/ajaffe/Amanda/libd_stem_timecourse/IRFinder/hg38"
sampIDs = as.vector(pd$SampleID)

#### intron map ####

# intronMap = read.table(file.path(idir, sampIDs[1], "IRFinder-IR-dir.txt"),header = TRUE)[,c(1:3,6)]
# intronMap$Length = intronMap$End-intronMap$Start+1
# intronMap$intron_id = paste0("chr",intronMap$Chr,":",intronMap$Start,"-",intronMap$End,"(",intronMap$Direction,")")
# intronMap$Chr = paste0("chr",intronMap$Chr)

# ## add gene info according to coords
# library('BSgenome.Hsapiens.UCSC.hg38')
# gtf = import(con="/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf", format="gtf")
# genes = as.data.frame(gtf[which(gtf$type=="gene" & gtf$gene_status=="KNOWN")])[,1:13]
# genes$strand = droplevels(genes$strand)  ## drop "*" as factor level
# rm(gtf)
# intronMap$Gene = bpmapply(function(ch,s,e,d) 
		# paste(genes[which(genes$seqnames==ch & genes$start<s & genes$end>e & genes$strand==d),"gene_name"],collapse=","), 
		# intronMap$Chr, intronMap$Start, intronMap$End, intronMap$Direction, 
		# BPPARAM = MulticoreParam(8))

# ## add Entrez ID
# library('biomaRt')
# ensembl = useMart("ENSEMBL_MART_ENSEMBL",  
	# dataset="hsapiens_gene_ensembl", host="jul2016.archive.ensembl.org")
# sym = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene"), 
		# values=geneMap$ensemblID, mart=ensembl)
# intronMap$EntrezID = sym$entrezgene[match(intronMap$Gene, sym$hgnc_symbol)]

		
# ## other results		
# coverage = bplapply(sampIDs, function(x) {
	# read.table(file.path(idir, x, "IRFinder-IR-dir.txt"),header = TRUE)$Coverage }, 
	# BPPARAM = MulticoreParam(8))
# coverage = do.call(cbind,coverage)

# intronDepth = bplapply(sampIDs, function(x) {
	# read.table(file.path(idir, x, "IRFinder-IR-dir.txt"),header = TRUE)$IntronDepth }, 
	# BPPARAM = MulticoreParam(8))
# intronDepth = do.call(cbind,intronDepth)

# IRratio = bplapply(sampIDs, function(x) {
	# read.table(file.path(idir, x, "IRFinder-IR-dir.txt"),header = TRUE)$IRratio }, 
	# BPPARAM = MulticoreParam(8))
# IRratio = do.call(cbind,IRratio)

# rownames(IRratio) = intronMap$intron_id
# colnames(IRratio) = pd$RNA_NO

# save(intronMap, coverage, intronDepth, IRratio, file="intron_data_tables.rda")

load("intron_data_tables.rda")



####################################
# run modeling
library(lme4)
library(lmerTest)
library(doMC)
registerDoMC(cores=12)

pddf = as.data.frame(pd)

## Remove introns with mostly 0 or mostly 1 ratios (less than 3 different)
keepInd = which(rowSums(as.data.frame(IRratio!=0)) >= 3 & rowSums(as.data.frame(IRratio!=1)) <= 103 )
IRratio = IRratio[keepInd,]
intronMapIR = intronMap[keepInd,]
intronDepth = intronDepth[keepInd,]
coverage = coverage[keepInd,]

# ## using factor
# options(warn=-1)
# lmer_list = foreach(i = 1:nrow(IRratio)) %dopar% {
	# if(i %% 5000 == 0) cat(".")
	# fit = lm(IRratio[i,] ~ COND + LINE + totalAssignedGene, data=pddf)
	# list(ANOVA = anova(fit), COEF=summary(fit)$coef)
# }
# names(lmer_list) = rownames(IRratio)

# # ##check output
# # which(sapply(lmer_list,function(x) ncol(x$COEF)) != 4)
# # which(sapply(lmer_list,function(x) ncol(x$ANOVA)) != 5)

# save(lmer_list, file="lm_IRratio3.rda")


####################################
### analysis ###########

load("lm_IRratio3.rda")


#### anova
anovaList = lapply(lmer_list, function(x) x$ANOVA)
anovaList = mclapply(anovaList, function(x) as.matrix(x)[1,4:5],mc.cores=12)
anovaMat= do.call("rbind", anovaList)

outStats = as.data.frame(anovaMat)
colnames(outStats) = c("anovaF", "anovaPvalue")
outStats$anovaBonf = p.adjust(outStats$anovaPvalue,"bonf")
outStats$anovaFDR = p.adjust(outStats$anovaPvalue,"fdr")

## linear
slopeMat = t(sapply(lmer_list, 
	function(x) x$COEF[2:4,1]))
colnames(slopeMat) = paste0("log2FC_",levels(pd$COND)[-1])
tMat = t(sapply(lmer_list, 
	function(x) x$COEF[2:4,3]))
colnames(tMat) = paste0("t_",levels(pd$COND)[-1])
pvalMat = t(sapply(lmer_list, 
	function(x) x$COEF[2:4,4]))
colnames(pvalMat) = paste0("pvalue_",levels(pd$COND)[-1])
qvalMat= matrix(, nrow = nrow(pvalMat), ncol = 3)
qvalMat[,1:3] = p.adjust(pvalMat[,1:3],method="fdr") 
colnames(qvalMat) = paste0("qvalue_",levels(pd$COND)[-1])

outStats = cbind(outStats, slopeMat,tMat,pvalMat,qvalMat)
rownames(outStats) = rownames(intronMapIR)


## add annotation
out = cbind(intronMapIR, outStats)

sum(out$anovaFDR < 0.05)
# 6285


#### csv of top results
ffSig = out[order(out$anovaPvalue),]
ffSig = ffSig[ffSig$anovaFDR < 0.05,]
ffSig$AD_meanIR = signif(rowMeans(IRratio[ffSig$intron_id,which(pd$CONDITION=="ACC_DORSAL(2)")] ),2)
ffSig$NPC_meanIR = signif(rowMeans(IRratio[ffSig$intron_id,which(pd$CONDITION=="NPC")] ),2)
ffSig$ROSE_meanIR = signif(rowMeans(IRratio[ffSig$intron_id,which(pd$CONDITION=="ROSETTE")] ),2)
ffSig$NEU_meanIR = signif(rowMeans(IRratio[ffSig$intron_id,which(pd$CONDITION=="NEURONS_PLUS_ASTROS")] ),2)
write.csv(ffSig, file="intronRatios_top_effects.csv", quote=FALSE)



####################################
### plot ###########

# Colors by Donor
pal = c(colorRampPalette(brewer.pal(10,"Paired")[1:2])(4), 
		colorRampPalette(brewer.pal(10,"Paired")[3:4])(3), 
		brewer.pal(10,"Paired")[9:10], brewer.pal(10,"Paired")[5:8])
pd$lineCol = pal[as.numeric(factor(pd$LINE))]


sigOrderMat = as.data.frame(apply(out[,10:11], 2, 
                                  function(x) order(x)[1:100]))

ooL = sigOrderMat$"anovaPvalue"
pdf("intronRatios_top_effects_cov.pdf",h=12,w=16)
par(mfrow=c(2,2),mar=c(5,6,5,2),cex.axis=1.5,cex.lab=2,cex.main=2)
for(i in ooL) {
 plot(IRratio[i,] ~ jitter(as.numeric(pd$DAY),.5), xaxt="n",
	   pch = 21, bg=pd$lineCol,
       cex=2,xlab="Day",
       ylab="IR Ratio", ylim=c(0,1.15),
       main = paste0(intronMapIR$intron_id[i], "\n", intronMapIR$Gene[i]) )
 axis(1, at=1:1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
 abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
 legend("top", paste0("p=",signif(out$"anovaPvalue"[i],3)), bg="white")
 if (i==ooL[1]) { legend("topright", levels(factor(pd$LINE)),
	pch = 15, col = pal, cex=1) }
	
 plot(intronDepth[i,] ~ jitter(as.numeric(pd$DAY),.5), xaxt="n",
	   pch = 21, bg=pd$lineCol,
       cex=2,xlab="Day", 
       ylab="Depth")
 axis(1, at=1:1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
 abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
 
 plot(1, type="n", axes=F, xlab="", ylab="")
 
  plot(coverage[i,] ~ jitter(as.numeric(pd$DAY),.5), xaxt="n",
	   pch = 21, bg=pd$lineCol,
       cex=2,xlab="Day", ylim=c(0,1),
       ylab="Coverage")
 axis(1, at=1:1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
 abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
}
dev.off()



### mean by day ###########

IR2 = colMeans(IRratio)

## across all introns
pdf("intronRatios_Overall.pdf",h=6,w=8)
par(mar=c(5,6,5,2),cex.axis=1.8,cex.lab=2,cex.main=2)
 boxplot(IR2 ~ pd$DAY,
          ylim = c(.1,.4),
          outline=FALSE, ylab="Mean IR Ratio", xlab="Day")
  points(IR2 ~ jitter(as.numeric(factor(pd$DAY)), 
		  amount=0.15), cex=1.5,
		  pch=21, bg = pd$lineCol)
dev.off()




################
# associations #
################

# gene set associations
library(clusterProfiler)
library(org.Hs.eg.db)

# split genes into up/down regulated, dropping grey
moduleGeneList = out$EntrezID[out$anovaBonf < 0.05]
moduleGeneList = moduleGeneList[!is.na(moduleGeneList)]

moduleGeneList = list(moduleGeneList)
names(moduleGeneList) = c("IRratio")

## set universe of expressed genes
geneUniverse = as.character(out$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

############################## 
## run enrichment analysis ###
##############################

## adjusted
goBP <- compareCluster(moduleGeneList, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Hs.eg.db,
                ont = "BP", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
				readable= TRUE)
goMF <- compareCluster(moduleGeneList, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Hs.eg.db,
                ont = "MF", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
				readable= TRUE)
goCC <- compareCluster(moduleGeneList, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Hs.eg.db,
                ont = "CC", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
				readable= TRUE)
kegg <- compareCluster(moduleGeneList, fun = "enrichKEGG",
                universe = geneUniverse,  pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1)
save(goBP,goMF,goCC,kegg, file="compareClusters_introns_bonf.rda")

load("compareClusters_introns.rda")

# ### combine
# geneSetStats_Adj = rbind(as.data.frame(goBP_Adj), as.data.frame(goMF_Adj),
	# as.data.frame(goCC_Adj), as.data.frame(kegg_Adj))
# geneSetStats_Adj$Ontology = rep(c("BP", "MF","CC", "KEGG"), 
	# times = c(nrow(goBP_Adj), nrow(goMF_Adj), nrow(goCC_Adj), nrow(kegg_Adj)))
# geneSetStats_Adj$Cluster = paste0("Adj_", geneSetStats_Adj$Cluster)

pdf("intron_ratio_enrichments_stemcell_bonf.pdf",h=5,w=7)
dotplot(goBP, includeAll="FALSE")
dotplot(goMF, includeAll="FALSE")
dotplot(goCC, includeAll="FALSE")
dotplot(kegg, includeAll="FALSE")
dev.off()



################
# separate up/down #
################

# split genes into condition lists
up = out$EntrezID[out$anovaFDR < 0.001 & out$log2FC_NEURONS_PLUS_ASTROS>0]
down = out$EntrezID[out$anovaFDR < 0.001 & out$log2FC_NEURONS_PLUS_ASTROS<0]
moduleGeneList = list(up, down)
moduleGeneList = lapply(moduleGeneList, function(x) x[!is.na(x)])
names(moduleGeneList) = c("IRratio-Increasing","IRratio-Decreasing")

## set universe of expressed genes
geneUniverse = as.character(out$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

############################## 
## run enrichment analysis ###
##############################

## adjusted
goBP <- compareCluster(moduleGeneList, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Hs.eg.db,
                ont = "BP", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
				readable= TRUE)
goMF <- compareCluster(moduleGeneList, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Hs.eg.db,
                ont = "MF", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
				readable= TRUE)
goCC <- compareCluster(moduleGeneList, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Hs.eg.db,
                ont = "CC", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
				readable= TRUE)
kegg <- compareCluster(moduleGeneList, fun = "enrichKEGG",
                universe = geneUniverse,  pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1)
save(goBP,goMF,goCC,kegg, file="compareClusters_introns_updown_fdr001.rda")

load("compareClusters_introns_updown_fdr001.rda")

# ### combine
# geneSetStats_Adj = rbind(as.data.frame(goBP_Adj), as.data.frame(goMF_Adj),
	# as.data.frame(goCC_Adj), as.data.frame(kegg_Adj))
# geneSetStats_Adj$Ontology = rep(c("BP", "MF","CC", "KEGG"), 
	# times = c(nrow(goBP_Adj), nrow(goMF_Adj), nrow(goCC_Adj), nrow(kegg_Adj)))
# geneSetStats_Adj$Cluster = paste0("Adj_", geneSetStats_Adj$Cluster)

pdf("intron_ratio_enrichments_stemcell_updown_fdr001.pdf",h=5,w=8)
dotplot(goBP, includeAll="FALSE")
dotplot(goMF, includeAll="FALSE")
dotplot(goCC, includeAll="FALSE")
dotplot(kegg, includeAll="FALSE")
dev.off()






############################## 
## check top results exons ###
##############################

top = out[sigOrderMat$anovaPvalue,]

######################
## load data
MAINDIR = '/dcl01/lieber/ajaffe/lab/libd_stem_timecourse'
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseExon_n146.rda"))

exonCounts = assays(rse_exon)$counts
exonMap = as.data.frame(rowRanges(rse_exon))
exonMap$strand = droplevels(exonMap$strand)
exonMap$eID = paste0(exonMap$seqnames,":", exonMap$start,"-", exonMap$end,"(",exonMap$strand,")" )

## recount getRPKM version
getRPKM <- function(rse, length_var = 'bp_length', mapped_var = NULL) {
    mapped <- if(!is.null(mapped_var)) colData(rse)[, mapped_var] else colSums(assays(rse)$counts)      
    bg <- matrix(mapped, ncol = ncol(rse), nrow = nrow(rse), byrow = TRUE)   
    len <- if(!is.null(length_var)) rowData(rse)[, length_var] else width(rowRanges(rse))   
    wid <- matrix(len, nrow = nrow(rse), ncol = ncol(rse), byrow = FALSE)
    assays(rse)$counts / (wid/1000) / (bg/1e6)
}
exonRpkm = getRPKM(rse_exon, length_var="Length")

exonCounts = exonCounts[,pd$SampleID]
exonRpkm = exonRpkm[,pd$SampleID]
yExprs = log2(exonRpkm+1)

top$exonUpind = top$exonDownind = NA
for (i in 1:100) {
top$exonUpind[i] = paste(which(exonMap$end==top$Start[i] & exonMap$strand==top$Direction[i]),collapse=",")
top$exonDownind[i] = paste(which(exonMap$start==top$End[i]+1 & exonMap$strand==top$Direction[i]),collapse=",")
}



pdf("intronRatios_exons.pdf",h=5,w=16)
par(mfrow=c(1,3),mar=c(5,6,5,2),cex.axis=1.5,cex.lab=2,cex.main=2)
for(i in intersect(grep(",",top$exonUpind, invert=T), grep(",",top$exonDownind, invert=T))) {
 plot(yExprs[as.numeric(top$exonUpind[i]),] ~ jitter(as.numeric(pd$DAY),.5), xaxt="n",
	   pch = 21, bg=pd$lineCol,
       cex=2,xlab="Day",
       ylab="log2(Exprs+1)", 
       main = paste0("exon\n",exonMap$eID[as.numeric(top$exonUpind[i])]) )
 axis(1, at=1:1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
 abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
 
 plot(IRratio[top$intron_id[i],] ~ jitter(as.numeric(pd$DAY),.5), xaxt="n",
	   pch = 21, bg=pd$lineCol,
       cex=2,xlab="Day",
       ylab="IR Ratio", ylim=c(0,1.15),
       main = paste0(top$intron_id[i], "\n", top$Gene[i]) )
 axis(1, at=1:1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
 abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
 if (i==1) { legend("topright", levels(factor(pd$LINE)),
	pch = 15, col = pal, cex=1) }

  plot(yExprs[as.numeric(top$exonDownind[i]),] ~ jitter(as.numeric(pd$DAY),.5), xaxt="n",
	   pch = 21, bg=pd$lineCol,
       cex=2,xlab="Day",
       ylab="log2(Exprs+1)", 
       main = paste0("exon\n",exonMap$eID[as.numeric(top$exonDownind[i])]) )
 axis(1, at=1:1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
 abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
}
dev.off()






###########################################
###########################################
###########################################
###########################################
###########################################
###########################################
### Compare IR ratios of neurons alone, neurons plus astros, and npc/rosettes

### load packages
library(jaffelab)
library(GenomicRanges)
library(rtracklayer)
library(SummarizedExperiment)
library(RColorBrewer)
library(rafalib)
library(edgeR)
library(limma)
library(BiocParallel)
library(recount)

######################
## load data
MAINDIR = '/dcl01/lieber/ajaffe/lab/libd_stem_timecourse'
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n146.rda"))

geneCounts = assays(rse_gene)$counts
geneMap = as.data.frame(rowRanges(rse_gene))

geneRpkm = recount::getRPKM(rse_gene, length_var="Length")


pd = colData(rse_gene)
pd$DAY = ordered(pd$DAY, levels=sort(as.numeric(levels(as.factor(pd$DAY)))) )
pd$DIV = as.numeric(as.character(pd$DAY))
pd$EXPERIMENT[pd$EXPERIMENT=="TIMECOURSE1"] = "TIMECOURSE"
pd$LINE = gsub('-','_',pd$LINE)

keepInd = which(pd$CONDITION %in% c("NEURONS_ALONE","NEURONS_PLUS_ASTROS","NPC","ROSETTE"))
pd = pd[keepInd,]
geneRpkm = geneRpkm[,keepInd]
geneCounts = geneCounts[,keepInd]   # n=106

######################
######################

intronDir = "/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/Counts/intron"
intronFn = file.path(intronDir, paste0(pd$SAMPLE_ID, '_Gencode.v25.hg38_Introns.counts.summary'))

# number of reads assigned to introns
intronStatList = lapply(intronFn, read.delim,row.names=1)
intronStats = do.call("cbind", intronStatList)
colnames(intronStats) = pd$SAMPLE_ID

pd$assignedIntron = as.numeric(intronStats[1,])
pd$assignedIntronRate = as.numeric(intronStats[1,] / colSums(intronStats))
pd$assignedIntronRate2 = as.numeric(intronStats[1,] / (intronStats[1,]+intronStats[2,]+intronStats[4,]))

## p-value of difference between npc/rosette and neuron_plus intron rates
pd$COND = pd$CONDITION
pd$COND[which(pd$COND %in% c("NPC","ROSETTE"))] = "NPCROSETTE"

fit = lm(assignedIntronRate2 ~ COND + LINE, data=pd)
summary(fit)
summary(fit)$coef[2,4] ## neurons_alone vs. neurons_plus_astros
summary(fit)$coef[3,4] ## neurons_alone vs. npc/rosette





