###

library(jaffelab)
library(limma)
library(GenomicRanges)
library(sva)
library(GOstats)
library(WGCNA)
library(SummarizedExperiment)
library(RColorBrewer)

## WGCNA options
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads=8)

######################
## load data
MAINDIR = '/dcl01/lieber/ajaffe/lab/libd_stem_timecourse'
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n146.rda"))

## recount getRPKM version
getRPKM <- function(rse, length_var = 'bp_length', mapped_var = NULL) {
    mapped <- if(!is.null(mapped_var)) colData(rse)[, mapped_var] else colSums(assays(rse)$counts)      
    bg <- matrix(mapped, ncol = ncol(rse), nrow = nrow(rse), byrow = TRUE)   
    len <- if(!is.null(length_var)) rowData(rse)[, length_var] else width(rowRanges(rse))   
    wid <- matrix(len, nrow = nrow(rse), ncol = ncol(rse), byrow = FALSE)
    assays(rse)$counts / (wid/1000) / (bg/1e6)
}
geneRpkm = getRPKM(rse_gene, length_var="Length")

geneCounts = assays(rse_gene)$counts
geneMap = rowData(rse_gene)
geneMap = as.data.frame(rowRanges(rse_gene))
# geneRpkm = getRPKM(rse_gene)

pd = colData(rse_gene)
pd$EXPERIMENT[pd$EXPERIMENT=="TIMECOURSE1"] = "TIMECOURSE"
pd$LINE = gsub('-','_',pd$LINE)
pd$DAY = ordered(pd$DAY, levels=sort(as.numeric(levels(as.factor(pd$DAY)))) )


# Colors by Donor
pal = c(colorRampPalette(brewer.pal(10,"Paired")[1:2])(4), 
		colorRampPalette(brewer.pal(10,"Paired")[3:4])(3), 
		brewer.pal(10,"Paired")[9:10], brewer.pal(10,"Paired")[5:8])
pd$lineCol = pal[as.numeric(factor(pd$LINE))]


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

######################

## drop low expression
geneIndex = which(rowMeans(geneRpkm) > 0.1)
geneRpkm = geneRpkm[geneIndex,]
geneCounts = geneCounts[geneIndex,]
geneMap = geneMap[geneIndex,]
## log
geneExprs = log2(geneRpkm+1) # transform
# ## or DESeq2 transform
# library(DESeq2)
# geneExprs <- DESeq2::varianceStabilizingTransformation(geneCounts)

mod = model.matrix(~COND + LINE + totalAssignedGene, data=pd)

## normalize
geneExprsAdj = cleaningY(geneExprs, mod, P=16)

################
## thresholding
# powers = c(c(1:10), seq(from = 12, to=20, by=2))
# sftAdj = pickSoftThreshold(t(geneExprsAdj), networkType = "signed",
	# powerVector = powers, verbose = 5)
# save(sftAdj, file="rdas/wgcna_soft_threshold_stemcell_signed.rda")
load("rdas/wgcna_soft_threshold_stemcell_signed.rda")

## threshold
sftAdj$powerEstimate
## = 6

# pdf("wgcna_power_recount.pdf")
# plot(sftAdj$fitIndices[,1], sftAdj$fitIndices[,2], 
	# type = "b",ylim=c(0,1))
# plot(sftAdj$fitIndices[,1], sftAdj$fitIndices[,5], 
	# type = "b",ylim=c(0,5000))
# dev.off()

# pdf('wgcna_power_recount_signed.pdf')
# par(mfrow = c(1,2)); cex1 = 0.9
# # Scale-free topology fit index as a function of the soft-thresholding power
# plot(sftAdj$fitIndices[,1], -sign(sftAdj$fitIndices[,3])*sftAdj$fitIndices[,2],
     # xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     # main = paste("Scale independence"))
# text(sftAdj$fitIndices[,1], -sign(sftAdj$fitIndices[,3])*sftAdj$fitIndices[,2],
     # labels=powers,cex=cex1,col="red")
# # this line corresponds to using an R^2 cut-off of h
# abline(h=0.90,col="red")

# plot(sftAdj$fitIndices[,1], sftAdj$fitIndices[,5],
     # xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     # main = paste("Mean connectivity"))
# text(sftAdj$fitIndices[,1], sftAdj$fitIndices[,5], labels=powers, cex=cex1,col="red")
# dev.off()
	
print("Starting clustering")	
##############	
## clustering
# netAdj = blockwiseModules(t(geneExprsAdj), power = sftAdj$powerEstimate,
	# TOMType = "signed", networkType = "signed",
	# minModuleSize = 30,
	# reassignThreshold = 0, mergeCutHeight = 0.25,
	# numericLabels = TRUE, pamRespectsDendro = FALSE,
	# saveTOMs = TRUE, 
	# saveTOMFileBase = "/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/tc_analysis/wgcna/LIBD_Adj_signed",
	# verbose = 3, maxBlockSize = 10000)
# save(netAdj, file = "rdas/stemcell_modules_signed.rda")
load("rdas/stemcell_modules_signed.rda")

table(netAdj$colors)
	# signed
  # 0    1    2    3    4    5    6    7    8    9   10   11
# 3284 8873 5942 2027 1471 1394 1213  592  256  249   94   71



 
################
####   pca  ####
################

rownames(geneMap) = geneMap$gencodeID

moduleInds = split(1:nrow(geneRpkm), netAdj$colors)
names(moduleInds) = paste0("Cluster_", names(moduleInds))
moduleInds$"Cluster_0" = NULL

geneRpkmModules = lapply(moduleInds, function(x) geneRpkm[x,])
pcaModules = lapply(geneRpkmModules, function(x) prcomp(t(log2(x+1))) )


#
# Plot 200 with largest rotation for each cluster
for (k in 1:11) {

pca_k = abs(pcaModules[[k]]$rotation[,1:2])
sigOrderMat = data.frame(PC1 = order(pca_k[,1],decreasing=TRUE)[1:200])
sigOrderMat$gene = rownames(pca_k)[sigOrderMat$PC1]
sigOrderMat$Symbol = geneMap[sigOrderMat$gene,"Symbol"]

maxi = min(100, table(netAdj$colors)[k+1] )

## Clusters
pdf(paste0("cluster",k,"_topRotation.pdf"),h=6,w=8)
par(mar=c(5,6,5,2),cex.axis=1.5,cex.lab=2,cex.main=2)
for (i in 1:maxi) {
gene = sigOrderMat$gene[i]
symbol = sigOrderMat$Symbol[i]
plot(geneExprs[gene,] ~ jitter(as.numeric(pd$DAY),.5), xaxt="n",
    pch = 21, bg=pd$lineCol,
        cex=2,xlab="Day",
        ylab="log2(Exprs+1)",
		main=paste0(symbol,"\n",gene) )
  axis(1, at=1:1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
  abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
}
dev.off()

}

## Gene list tables by cluster
for (k in 1:6) {
gene_k = rownames(geneRpkmModules[[k]])
genelist = data.frame(gencodeID = gene_k)
genelist$wgcnaCluster = k
genelist$Chr = geneMap[genelist$gencodeID,"seqnames"]
genelist$Start = geneMap[genelist$gencodeID,"start"]
genelist$End = geneMap[genelist$gencodeID,"end"]
genelist$Strand = geneMap[genelist$gencodeID,"strand"]
genelist$gene_type = geneMap[genelist$gencodeID,"gene_type"]
genelist$Symbol = geneMap[genelist$gencodeID,"Symbol"]

nam = paste0("cluster",k)
assign(nam, genelist)
write.csv(genelist, file=paste0("cluster",k,"_genelist.csv") )
}



################
# associations #
################

# gene set associations
library(clusterProfiler)
library(org.Hs.eg.db)

# split genes into modules, dropping grey
moduleGeneList_adj = split(geneMap$EntrezID, netAdj$colors)
moduleGeneList_adj = lapply(moduleGeneList_adj, function(x) x[!is.na(x)])
moduleGeneList_adj = moduleGeneList_adj[-1]

## set universe of expressed genes
geneUniverse = as.character(geneMap$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

############################## 
## run enrichment analysis ###
##############################

# ## and adjusted
# goBP_Adj <- compareCluster(moduleGeneList_adj, fun = "enrichGO",
                # universe = geneUniverse, OrgDb = org.Hs.eg.db,
                # ont = "BP", pAdjustMethod = "BH",
                # pvalueCutoff  = .1, qvalueCutoff  = 1,
				# readable= TRUE)
# goMF_Adj <- compareCluster(moduleGeneList_adj, fun = "enrichGO",
                # universe = geneUniverse, OrgDb = org.Hs.eg.db,
                # ont = "MF", pAdjustMethod = "BH",
                # pvalueCutoff  = .1, qvalueCutoff  = 1,
				# readable= TRUE)
# goCC_Adj <- compareCluster(moduleGeneList_adj, fun = "enrichGO",
                # universe = geneUniverse, OrgDb = org.Hs.eg.db,
                # ont = "CC", pAdjustMethod = "BH",
                # pvalueCutoff  = .1, qvalueCutoff  = 1,
				# readable= TRUE)
# kegg_Adj <- compareCluster(moduleGeneList_adj, fun = "enrichKEGG",
                # universe = geneUniverse,  pAdjustMethod = "BH",
                # pvalueCutoff  = .1, qvalueCutoff  = 1)
# save(goBP_Adj, goMF_Adj, goCC_Adj, kegg_Adj, file="wgcna_compareCluster_signed.rda")
load("wgcna_compareCluster_signed.rda")
write.csv(as.data.frame(goBP_Adj), file="enrich_BP.csv")
write.csv(as.data.frame(goMF_Adj), file="enrich_MF.csv")
write.csv(as.data.frame(goCC_Adj), file="enrich_CC.csv")
write.csv(as.data.frame(kegg_Adj), file="enrich_kegg.csv")
				
pdf("wgcna_cluster_enrichments_stemcell_signed.pdf",h=6,w=10)
dotplot(goBP_Adj, includeAll="TRUE")
dotplot(goMF_Adj, includeAll="TRUE")
dotplot(goCC_Adj, includeAll="TRUE")
dotplot(kegg_Adj, includeAll="TRUE")
dev.off()				
		

		
################
# plot pattern #
################

datME = moduleEigengenes(t(geneExprsAdj), netAdj$colors)$eigengenes
datME = datME[,-1] ## remove unassigned

datMEmean = list()
for (i in 1:11) {
datMEmean[[i]] = unlist(lapply(split(datME[,i], pd$DAY), mean))
}

pdf("wgcna_patterns_eigengenes_signed.pdf",h=8,w=12)
par(mfrow=c(2,3), mar=c(5,5,5,2),cex.axis=1.3,cex.lab=1.8,cex.main=2)
for (i in 1:11) {
 plot(1:9, datMEmean[[i]], cex=2.5, type="b", xaxt="n", ylim=c(-.22,.22), 
		  pch=16, col="darkblue",
          ylab="Eigengene", xlab="Day", main=paste0("Cluster ",i,"  (", table(netAdj$colors)[i+1],")") )
 axis(1, at=1:length(levels(pd$DAY)), labels = levels(pd$DAY))
 abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
}
dev.off()


################
# all on one plot #

## plot legend next to plot
m <- t(as.matrix(c(1,1,1,1,1,2)))
mycols <- colors()[c(551,556,136,526,420,625,601,87,50,81,280)]

## 133,101,544,385,230,401,88,12,563,597,287
#### all 11 clusters in same panel
pdf("figure2_wgcna_eigengenes.pdf",h=5,w=6)
layout(m)
par(mar=c(5,5,5,1),cex.axis=2,cex.lab=2,cex.main=2)

palette(mycols)

i=1
plot(1:9, datMEmean[[i]], cex=2.5, type="l", lwd=2, xaxt="n", ylim=c(-.22,.22), 
	  pch=16, col=i,
	  ylab="Eigengene", xlab="Day", main="WGCNA Modules" )
axis(1, at=1:length(levels(pd$DAY)), labels = levels(pd$DAY))
abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
for (i in 2:11) {
 points(1:9, datMEmean[[i]], cex=2.5, type="l", lwd=3, pch=16, col=i)
}
par(mar=c(0,0,5,1),cex.axis=2,cex.lab=2,cex.main=2)
plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
legend("topleft", paste0("Module ", 1:11),
	pch = 15, col = 1:11, cex=1.1, bty="n",pt.cex=2)

dev.off()



## Four panels
mycols <- colors()[c(81,30,618,556,134,551,628,601,50,87,421)]

pdf("figure2_wgcna_eigengenes_4_log.pdf",h=12,w=14)
par(mfrow=c(2,2),mar=c(5,5,1,1),cex.axis=2,cex.lab=2.5,cex.main=3)
palette(mycols)

i=1
plot(1:9, datMEmean[[i]], cex=1, type="b", lwd=1, xaxt="n", ylim=c(-.22,.22), 
	  pch=16, col=i,
	  ylab="Eigengene", xlab="", main="" )
axis(1, at=1:length(levels(pd$DAY)), labels = levels(pd$DAY), cex.axis=2.2)
abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
legend("bottomright", paste0("Module ", c(1,9,10)),
	pch = 15, col = c(1,9,10), cex=2, bty="n")
for (i in c(10,9,1)) {
 points(1:9, datMEmean[[i]], cex=log10(table(netAdj$colors)[i+1]), type="b", pch=16, col=i, lwd=log10(table(netAdj$colors)[i+1]) )
}

i=2
plot(1:9, datMEmean[[i]], cex=1, type="b", lwd=1, xaxt="n", ylim=c(-.22,.22), 
	  pch=16, col=i,
	  ylab="", xlab="", main="" )
axis(1, at=1:length(levels(pd$DAY)), labels = levels(pd$DAY), cex.axis=2.2)
abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
legend("bottomleft", paste0("Module ", c(2,3,8)),
	pch = 15, col = c(2,3,8), cex=2, bty="n")
for (i in c(8,3,2)) {
 points(1:9, datMEmean[[i]], cex=log10(table(netAdj$colors)[i+1]), type="b", pch=16, col=i, lwd=log10(table(netAdj$colors)[i+1]) )
}

i=4
plot(1:9, datMEmean[[i]], cex=1, type="b", lwd=1, xaxt="n", ylim=c(-.22,.22), 
	  pch=16, col=i,
	  ylab="Eigengene", xlab="Day", main="" )
axis(1, at=1:length(levels(pd$DAY)), labels = levels(pd$DAY), cex.axis=2.2)
abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
legend("bottomright", paste0("Module ", c(4,5,11)),
	pch = 15, col = c(4,5,11), cex=2, bty="n")
for (i in c(11,5,4)) {
 points(1:9, datMEmean[[i]], cex=log10(table(netAdj$colors)[i+1]), type="b", pch=16, col=i, lwd=log10(table(netAdj$colors)[i+1]) )
}

i=6
plot(1:9, datMEmean[[i]], cex=1, type="b", lwd=1, xaxt="n", ylim=c(-.22,.22), 
	  pch=16, col=i,
	  ylab="", xlab="Day", main="" )
axis(1, at=1:length(levels(pd$DAY)), labels = levels(pd$DAY), cex.axis=2.2)
abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
legend("bottomleft", paste0("Module ", c(6,7)),
	pch = 15, col = c(6,7), cex=2, bty="n")
for (i in c(7,6)) {
 points(1:9, datMEmean[[i]], cex=log10(table(netAdj$colors)[i+1]), type="b", pch=16, col=i, lwd=log10(table(netAdj$colors)[i+1]) )
}

dev.off()





# datMEmed = list()
# for (i in 1:9) {
# datMEmed[[i]] = unlist(lapply(split(datME[,i], pd$DAY), median))
# }
# pdf("wgcna_patterns.pdf",h=8,w=18)
# par(mfrow=c(2,5), mar=c(5,5,5,2),cex.axis=1.3,cex.lab=1.8,cex.main=2)
# for (i in 1:9) {
 # plot(1:9, datMEmed[[i]], cex=2.5, type="b", xaxt="n", ylim=c(-.22,.22), 
		  # pch=16, col="darkblue",
          # ylab="Eigengene", xlab="Day", main=paste0("Cluster ",i,"  (", table(netAdj$colors)[i+1],")") )
 # axis(1, at=1:length(levels(pd$DAY)), labels = levels(pd$DAY))
 # abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
# }
# dev.off()



# ## Clusters
# for (clus in 1:9) {
 # pdf(paste0("cluster",clus,".pdf"),h=6,w=8)
 # par(mar=c(5,6,5,2),cex.axis=1.5,cex.lab=2,cex.main=2)
 # for(i in which(netAdj$colors==clus)[1:min(200,table(netAdj$colors)[clus+1])] ) {
  # plot(geneRpkm[i,] ~ jitter(as.numeric(pd$DAY),.5), xaxt="n",
    # pch = 21, bg=pd$lineCol,
        # cex=2,xlab="Day",
        # ylab="log2(Exprs+1)",
		# main=geneMap$Symbol[i] )
  # axis(1, at=1:1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
  # abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
# }
# dev.off()
# }



# ################
# #  summarize   #
# ################

# ### combine
# geneSetStats_Adj = rbind(as.data.frame(goBP_Adj), as.data.frame(goMF_Adj),
	# as.data.frame(goCC_Adj), as.data.frame(kegg_Adj))
# geneSetStats_Adj$Ontology = rep(c("BP", "MF","CC", "KEGG"), 
	# times = c(nrow(goBP_Adj), nrow(goMF_Adj), nrow(goCC_Adj), nrow(kegg_Adj)))
# geneSetStats_Adj$Cluster = paste0("Adj_", geneSetStats_Adj$Cluster)
# geneSetStats = geneSetStats_Adj	

# ## reduce
# geneSetPvals = matrix(NA, nr = length(unique(geneSetStats$ID)),
	# nc = length(unique(geneSetStats$Cluster)), dimnames = list(
	# unique(geneSetStats$ID), unique(geneSetStats$Cluster)))
# geneSetList = split(geneSetStats, geneSetStats$Cluster)[colnames(geneSetPvals)]

# for(i in seq(along=geneSetList)) {
	# xx = geneSetList[[i]]
	# geneSetPvals[,i] = xx$pvalue[match(rownames(geneSetPvals), xx$ID)]
# }

# pvals = as.numeric(geneSetPvals)
# pvals[is.na(pvals)] = 1
# geneSetQvals = p.adjust(pvals, "fdr")
# geneSetQvals = matrix(geneSetQvals, nr = length(unique(geneSetStats$ID)),
	# nc = length(unique(geneSetStats$Cluster)), dimnames = list(
	# unique(geneSetStats$ID), unique(geneSetStats$Cluster)))
	
# ## filter 
# geneSetSig = geneSetQvals[rowSums(geneSetQvals < 0.05, na.rm=TRUE) > 0,]
# colSums(geneSetSig < 0.05, na.rm=TRUE)
# geneSetSig = as.data.frame(geneSetSig)
# geneSetSig$Description = geneSetStats$Description[
	# match(rownames(geneSetSig),geneSetStats$ID)]

# ## summarize
# geneSetSig$Summarized = apply(geneSetSig < 0.05, 1, function(x) 
	# paste(names(x[which(x)]),collapse=";"))
# save(geneSetSig, geneSetStats, file = "rdas/geneSetEnrichment_stemcell_WGCNA.rda")
# write.csv(geneSetSig, file = "geneSetEnrichment_stemcell_WGCNA.csv")


# #####################
# # Dx association ####
# #####################
# coefAdj = t(apply(netAdj$MEs, 2, function(y) 
	# summary(lm(y~ pd$Dx + pd$Age))$coef[2,]))
# coefAdjPosthoc = t(apply(netAdj$MEs, 2, function(y) 
	# summary(lm(y~ cbind(mod[,1:3], qSVs)- 1))$coef[2,]))
# rownames(coefAdj) =rownames(coefAdjPosthoc) = gsub("ME", "Adj_", rownames(coefAdj))
# coefAdj = coefAdj[paste0("Adj_", 0:(nrow(coefAdj)-1)),]
# coefAdjPosthoc = coefAdjPosthoc[paste0("Adj_", 0:(nrow(coefAdjPosthoc)-1)),]

# ## plots
# plot(-log10(coefAdj[,4]), -log10(coefAdjPosthoc[,4]),
	# main = "-log10(p-values) for Dx", xlab = "Adj Only",
	# ylab = "Posthost ME qSVA Adj")

# ## qSVA adjusted
# coefQsva = t(apply(netQsva$MEs, 2, function(y) 
	# summary(lm(y~pd$Dx + pd$Age))$coef[2,]))
# rownames(coefQsva) = gsub("ME", "qSVA_", rownames(coefQsva))
# coefQsva = coefQsva[paste0("qSVA_", 0:(nrow(coefQsva)-1)),]

# ## ME association with qSVs in unadjusted
# anovaConfound = apply(netAdj$MEs, 2, function(y) 
	# anova(lm(y~ qSVs))$"Pr(>F)"[1])
# names(anovaConfound)  = gsub("ME", "Adj_", names(anovaConfound))
# anovaConfound = anovaConfound[paste0("Adj_", 0:(length(anovaConfound)-1))]

# coefConfoundList = lapply(as.data.frame(netAdj$MEs), function(y) 
	# summary(lm(y~ qSVs))$coef)
# names(coefConfoundList)  = gsub("ME", "Adj_", names(coefConfoundList))
# coefConfoundList = coefConfoundList[names(anovaConfound)]

# ####################
# #### compare #######
# ####################

# netAdj_LIBD = netAdj
# netQsva_LIBD = netQsva
# geneMap_LIBD = geneMap
# pd_LIBD = pd
# mod_LIBD = mod
# qSVs_LIBD = qSVs

# #################
# ## within LIBD ##
# #################

# # get overlap
# tt_LIBD = table(netQsva_LIBD$colors, netAdj_LIBD$colors,
	# dnn = c("qSVA", "Adj"))
# tt_LIBD_bg = tt_LIBD
# tt_qsva = table(netQsva_LIBD$colors) ; tt_adj= table(netAdj_LIBD$colors)
# for(i in 0:max(netQsva_LIBD$colors)) {
	# for(j in 0:max(netAdj_LIBD$colors)) {
		# tt_LIBD_bg[as.character(i),as.character(j)] = tt_qsva[as.character(i)] + tt_adj[as.character(j)]
# }}
# overlap_LIBD = tt_LIBD / (tt_LIBD_bg - tt_LIBD)

# round(overlap_LIBD[-1,-1],2)
# apply(overlap_LIBD, 1, max)

# # coverage
# coverage_LIBD_qSVA = apply(prop.table(tt_LIBD,1),1,function(x) max(x[-1]))
# coverage_LIBD_Adj = apply(prop.table(tt_LIBD,2),2,function(x) max(x[-1]))
# plot(apply(overlap_LIBD, 1, function(x) max(x[-1])), coverage_LIBD_qSVA, type="n",
	# xlab = "Overlap", ylab = "Coverage (qSVA)")
# text(apply(overlap_LIBD, 1, function(x) max(x[-1])), coverage_LIBD_qSVA, 0:32)
# plot(apply(overlap_LIBD, 1, function(x) max(x[-1])), coverage_LIBD_qSVA, type="n",
	# xlab = "Overlap", ylab = "Coverage (qSVA)")
# text(apply(overlap_LIBD, 1, function(x) max(x[-1])), coverage_LIBD_qSVA, 0:32)

# #################
# # across data ###
# #################

# load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/wgcna/CMC_modules_twoRuns.rda")
# netAdj_CMC = netAdj
# netQsva_CMC = netQsva
# geneMap_CMC = geneMap
# load("rdas/geneLevel_CMC_qSVsAndMod.rda")
# pd_CMC = pd
# mod_CMC = mod
# qSVs_CMC = qSVs

# # match up genes
# table(rownames(geneMap_LIBD) %in% rownames(geneMap_CMC))
# mmGene = match(rownames(geneMap_LIBD), rownames(geneMap_CMC))

# ## get overlaps
# tt_qSVA = table(netQsva_LIBD$colors[!is.na(mmGene)], 
	# netQsva_CMC$colors[mmGene[!is.na(mmGene)]],	dnn = c("LIBD", "CMC"))
# tt_Cross_bg = tt_qSVA
# for(i in rownames(tt_qSVA)) {
	# for(j in colnames(tt_qSVA)) {
		# tt_Cross_bg[i,j] = rowSums(tt_qSVA)[i] + colSums(tt_qSVA)[j]
# }}
# overlap_Cross = tt_qSVA / (tt_Cross_bg - tt_qSVA)

# round(overlap_Cross[-1,-1],2)
# apply(overlap_Cross, 1, max)

# ## coverage
# coverage_Cross_qSVA = apply(prop.table(tt_qSVA,1),1,function(x) max(x[-1]))

# tt_Adj = table(netAdj_LIBD$colors[!is.na(mmGene)], 
	# netAdj_CMC$colors[mmGene[!is.na(mmGene)]],	dnn = c("LIBD", "CMC"))
# coverage_Cross_Adj = apply(prop.table(tt_Adj,1),1,function(x) max(x[-1]))

# # check plots
# plot(apply(overlap_Cross, 1, function(x) max(x[-1])), coverage_Cross_qSVA, type="n",
	# xlab = "Overlap", ylab = "Coverage (qSVA)")
# text(apply(overlap_Cross, 1, function(x) max(x[-1])), coverage_Cross_qSVA, 0:32)

# ## more diagnostics
# plot(coverage_LIBD_qSVA, coverage_Cross_qSVA, type="n",
	# xlab = "Within LIBD", ylab = "Across CMC", 
	# main = "qSVA Coverage",ylim = c(0,1),xlim = c(0,1))
# text(coverage_LIBD_qSVA, coverage_Cross_qSVA, 0:32)
# abline(0,1)

# plot(coverage_LIBD_Adj, coverage_Cross_Adj, type="n",
	# xlab = "Within LIBD", ylab = "Across CMC", 
	# main = "Adj Coverage",ylim = c(0,1),xlim = c(0,1))
# text(coverage_LIBD_Adj, coverage_Cross_Adj, 0:21)
# abline(0,1)

# plot(apply(overlap_LIBD, 1, function(x) max(x[-1])), 
	# apply(overlap_Cross, 1, function(x) max(x[-1])), type="n",
	# xlab = "Within LIBD", ylab = "Across CMC", 
	# main = "qSVA Overlap",ylim = c(0,1),xlim = c(0,1))
# text(apply(overlap_LIBD, 1, function(x) max(x[-1])), 
	# apply(overlap_Cross, 1, function(x) max(x[-1])), 0:32)
# abline(0,1)


# signif(prop.table(tt_qSVA[-1,-1],1),2)
# apply(signif(prop.table(tt_qSVA[-1,-1],1),2),1,max)
# apply(signif(prop.table(tt_qSVA[-1,-1],1),2),1,which.max)

# tt_Adj = table(netAdj_LIBD$colors[!is.na(mmGene)], 
	# netAdj_CMC$colors[mmGene[!is.na(mmGene)]],	dnn = c("LIBD", "CMC"))
# tt_Adj[-1,-1]
# signif(prop.table(tt_Adj[-1,-1],1),2)
# apply(signif(prop.table(tt_Adj[-1,-1],1),2),1,max)
# apply(signif(prop.table(tt_Adj[-1,-1],1),2),1,which.max)

# ## qsva across data versus within data 
# plot(apply(prop.table(tt_qSVA[-1,-1],1),1,max),
	# apply(prop.table(tt_LIBD[-1,-1],1),1,max),
	# xlab = "Across Data", ylab = "Across Method",
	# ylim = c(0,1),xlim = c(0,1))
# abline(0,1)

# # adj across data versus within LIBD
# plot(apply(prop.table(tt_Adj[-1,-1],1),1,max),
	# apply(prop.table(tt_LIBD[-1,-1],2),2,max),
	# xlab = "Across Data", ylab = "Across Method",
	# ylim = c(0,1),xlim = c(0,1))
# abline(0,1)

