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
load(file.path(MAINDIR,"rn6_pipe/rse_gene_timecourse_rat_n36.Rdata"))

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
geneMap = as.data.frame(rowRanges(rse_gene))

pd = colData(rse_gene)
pd$DIV = as.factor(pd$DAY)
pd$DIV = ordered(pd$DAY, levels=sort(as.numeric(levels(as.factor(pd$DAY)))) )


###########################
#### filter expression ####
geneMap$meanExprs = rowMeans(geneRpkm)
gIndex = which(geneMap$meanExprs > 0.1)
geneRpkm = geneRpkm[gIndex,]
geneCounts = geneCounts[gIndex,]
geneMap = geneMap[gIndex,]


############################
## do voom analysis ########

## don't include rat day 21 and 42 in model
keepInd = which(pd$DAY > 42)
pdSub = pd[keepInd,]
geneCountsSub = geneCounts[,keepInd]
geneRpkmSub = geneRpkm[,keepInd]

##### fix model
mod = model.matrix(~SPECIES + DAY , data=pdSub)

##############
## gene 
dge = DGEList(counts = geneCountsSub, genes = as.data.frame(geneMap))
dge = calcNormFactors(dge)
v = voom(dge, mod, plot=FALSE)
fit = lmFit(v)
ebGene = ebayes(fit)	
ffGene = topTable(eBayes(fit), coef=2, n = nrow(dge), adjust.method="fdr")
ffGene$bonf.P.Val = p.adjust(ffGene$P.Value, "bonf")
ffGene = ffGene[rownames(ebGene$t),]
ffGene$FC = ifelse(ffGene$logFC<0, "DOWN", "UP") ##DOWN = rat alone is lower

sum(ffGene$adj.P.Val < 0.05)
## 1329
sum(ffGene$bonf.P.Val < 0.05)



############################
## plot top ffGene

sigOrderMat = as.data.frame(apply(ffGene[,14:15], 2, 
                                  function(x) order(x)[1:200]))
yExprs = as.matrix(log2(geneRpkm+1))

## top effects plot
ooL = sigOrderMat$"P.Value"
pdf("rat_top_effects.pdf",h=6,w=6)
par(mar=c(5,6,5,2),cex.axis=1.4,cex.lab=2,cex.main=2)
palette(brewer.pal(8,"PuBuGn"))
for(i in ooL) {
 plot(yExprs[i,] ~ jitter(pd$DAY,.6), xaxt="n",
	   pch = 21, bg=as.numeric(pd$SPECIES)*4,
       cex=2,
	   xlab="Day", ylab="log2(Exprs + 1)", 
       main = paste0(geneMap$Symbol[i], "\n", geneMap$Geneid[i]) )
 axis(1, at=unique(pd$DAY), labels = unique(pd$DAY))
 legend("top", paste0("p=",signif(ffGene$"P.Value"[i],3)), bg="white")
 if (i==ooL[1]) { legend("bottomleft", levels(factor(pd$SPECIES)),
	pch = 15, col=c(4,8), cex=1) }
}
dev.off()


## volcano plot
pdf("rat_astrocytes_volcano.pdf",h=6,w=6)
par(mar=c(5,6,5,2),cex.axis=1.5,cex.lab=1.4,cex.main=1.5)
plot(ffGene$logFC, -log10(ffGene$P.Value),
	 xlim=c(-10,10),
	 xlab="log2(FoldChange)", ylab="-log10(p-value)",
	 pch=21, col="black",	 
	 bg=ifelse(ffGene$bonf.P.Val<0.05, "orangered", NA) ) 
# abline(v=c(seq(-15,15, by=5)), lty=2, col="grey")

plot(ffGene$logFC, -log10(ffGene$P.Value),
	 xlim=c(-10,10),
	 xlab="log2(FoldChange)", ylab="-log10(p-value)",
	 pch=21, col="black",	 
	 bg=ifelse(ffGene$adj.P.Val<0.05, "orangered", NA) )
legend("topleft", "FDR < 0.05", pch=16, col="orangered", cex=.8)
# abline(v=c(seq(-15,15, by=5)), lty=2, col="grey")
dev.off()



## pca
pca1 = prcomp(t(log2(geneRpkmSub+1)))
pcaVars1 = getPcaVars(pca1)

pdf("pca_log2Rpkm_rat.pdf", h=6,w=6)
par(mar=c(5,6,4,2),cex.axis=1.7,cex.lab=2,cex.main=2)
palette(brewer.pal(8,"PuBuGn"))
## PC 1 and 2
plot(pca1$x, bg=as.numeric(pdSub$SPECIES)*4, cex=2, main="Gene PCs",
	pch=21,
	xlab=paste0("PC1: ", pcaVars1[1], "% Var Expl"),
    ylab=paste0("PC2: ", pcaVars1[2], "% Var Expl"))
legend("topleft", paste0(levels(pdSub$SPECIES)), col=c(4,8), 
       pch = 15, cex=.9)
## PC 3 and 4
plot(pca1$x[,3],pca1$x[,4], bg=as.numeric(pdSub$SPECIES)*4, cex=2, main="Gene PCs",
	pch=21,
	xlab=paste0("PC3: ", pcaVars1[3], "% Var Expl"),
    ylab=paste0("PC4: ", pcaVars1[4], "% Var Expl"))
legend("topleft", paste0(levels(pdSub$SPECIES)), col=c(4,8), 
       pch = 15, cex=.9)

dev.off()	





################
# associations #
################

# gene set associations
library(clusterProfiler)
# library(org.Hs.eg.db)
library(org.Rn.eg.db)

identical(ffGene$Geneid, geneMap$Geneid) ##check order

# split genes into up/down regulated, dropping grey
moduleGeneList = geneMap$EntrezID[ffGene$adj.P.Val < 0.05]
moduleGeneList = moduleGeneList[!is.na(moduleGeneList)]

# split genes into condition lists
up = geneMap$EntrezID[ffGene$adj.P.Val < 0.05 & ffGene$FC=="UP"]
down = geneMap$EntrezID[ffGene$adj.P.Val < 0.05 & ffGene$FC=="DOWN"]
moduleGeneList = list(up, down)
moduleGeneList = lapply(moduleGeneList, function(x) x[!is.na(x)])
names(moduleGeneList) = c("Up-regulated","Down-regulated")
names(moduleGeneList) = c("Astro > Neuron+Astro","Astro < Neuron+Astro")

## set universe of expressed genes
geneUniverse = as.character(geneMap$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

############################## 
## run enrichment analysis ###
##############################

## adjusted
goBP <- compareCluster(moduleGeneList, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Rn.eg.db,
                ont = "BP", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = .01,
				readable= TRUE)
goMF <- compareCluster(moduleGeneList, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Rn.eg.db,
                ont = "MF", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = .01,
				readable= TRUE)
goCC <- compareCluster(moduleGeneList, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Rn.eg.db,
                ont = "CC", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = .01,
				readable= TRUE)
kegg <- compareCluster(moduleGeneList, fun = "enrichKEGG", organism="rno",
                universe = geneUniverse,  pAdjustMethod = "BH", 
                pvalueCutoff  = 1, qvalueCutoff  = .01)


# ### combine
# geneSetStats_Adj = rbind(as.data.frame(goBP_Adj), as.data.frame(goMF_Adj),
	# as.data.frame(goCC_Adj), as.data.frame(kegg_Adj))
# geneSetStats_Adj$Ontology = rep(c("BP", "MF","CC", "KEGG"), 
	# times = c(nrow(goBP_Adj), nrow(goMF_Adj), nrow(goCC_Adj), nrow(kegg_Adj)))
# geneSetStats_Adj$Cluster = paste0("Adj_", geneSetStats_Adj$Cluster)

pdf("ratAstro_up_down_enrichments_stemcell.pdf",h=5,w=9)
dotplot(goBP, includeAll="TRUE")
dotplot(goMF, includeAll="TRUE")
dotplot(goCC, includeAll="TRUE")
dotplot(kegg, includeAll="TRUE")
dev.off()




