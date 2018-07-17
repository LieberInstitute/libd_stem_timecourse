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
pd$REP = as.numeric(pd$BIO_REP)
pd$REP[is.na(pd$REP)] = 1
pd$DAY = ordered(pd$DAY, levels=sort(as.numeric(levels(as.factor(pd$DAY)))) )
pd$DIV = as.numeric(as.character(pd$DAY))
pd$EXPERIMENT[pd$EXPERIMENT=="TIMECOURSE1"] = "TIMECOURSE"

# Colors by Donor
pal = c(colorRampPalette(brewer.pal(10,"Paired")[1:2])(4), 
		colorRampPalette(brewer.pal(10,"Paired")[3:4])(3), 
		brewer.pal(10,"Paired")[9:10], brewer.pal(10,"Paired")[5:8])
pd$lineCol = pal[as.numeric(factor(pd$LINE))]

## only use NERUONS (day 77?)
bioInd = which(pd$Class == "Naked genomes" & grepl("NEURON",pd$CONDITION) &
                 pd$LINE %in% unique(pd$LINE[pd$CONDITION=='NEURONS_ALONE']) &
				 pd$DAY == 77)
pd = pd[bioInd,]
geneRpkm = geneRpkm[,bioInd]
geneCounts = geneCounts[,bioInd]

pd$COND = pd$CONDITION
pd$COND = as.factor(pd$COND) #alone v. plus astrocyte

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
mod = model.matrix(~COND + factor(LINE) + totalAssignedGene, data=pd)
colnames(mod)[2] = "ON_ASTROS"

##############
## gene 
dge = DGEList(counts = geneCounts, genes = as.data.frame(geneMap))
dge = calcNormFactors(dge)
v = voom(dge, mod, plot=FALSE)
fit = lmFit(v)
ebGene = ebayes(fit)	
ffGene = topTable(eBayes(fit), coef=2, n = nrow(dge), adjust.method="fdr")
ffGene$bonf.P.Val = p.adjust(ffGene$P.Value, "bonf")
ffGene = ffGene[rownames(ebGene$t),]
ffGene$FC = ifelse(ffGene$logFC<0, "DOWN", "UP")

sum(ffGene$adj.P.Val < 0.05)
# 3214
sum(ffGene$adj.P.Val < 0.01)
sum(ffGene$adj.P.Val < 0.005)
sum(ffGene$adj.P.Val < 0.001)

sum(ffGene$bonf.P.Val < 0.05)
sum(ffGene$bonf.P.Val < 0.01)
sum(ffGene$bonf.P.Val < 0.005)
sum(ffGene$bonf.P.Val < 0.001)


#### csv of top results
ffSig = ffGene[order(ffGene$P.Value),]
ffSig = ffSig[ffSig$adj.P.Val < 0.05,-which(colnames(ffSig)=="gencodeTx")]
ffSig$meanExprsOFF = rowMeans(geneRpkm[ffSig$gencodeID,which(pd$CONDITION=="NEURONS_ALONE")] )
ffSig$meanExprsON = rowMeans(geneRpkm[ffSig$gencodeID,which(pd$CONDITION=="NEURONS_PLUS_ASTROS")] )
write.csv(ffSig, file="astrocytes_day77_top_effects.csv", quote=FALSE)


# ############################
# ## plot top ffGene


# sigOrderMat = as.data.frame(apply(ffGene[,14:15], 2, 
                                  # function(x) order(x)[1:200]))
# yExprs = as.matrix(log2(geneRpkm+1))

# ## top effects plot
# ooL = sigOrderMat$"P.Value"
# pdf("astrocytes_day77_top_effects.pdf",h=6,w=6)
# par(mar=c(5,6,5,2),cex.axis=1.9,cex.lab=2,cex.main=2)
# for(i in ooL) {
 # plot(yExprs[i,] ~ jitter(as.numeric(pd$COND),.6) , xaxt="n",
	   # pch = ifelse(pd$DX=="CNT",21,23), bg=pd$lineCol,
       # cex=2, xlim=c(.75,2.25),
	   # xlab="Day 77", ylab="log2(Exprs + 1)", 
       # main = paste0(geneMap$Symbol[i], "\n", geneMap$gencodeID[i]) )
 # axis(1, at=1:2, labels = c("Off Astro","On Astro"))
 # legend("top", paste0("p=",signif(ffGene$"P.Value"[i],3)), bg="white")
 # if (i==ooL[1]) { legend("topleft", levels(factor(pd$LINE)),
	# pch = 15, col = levels(as.factor(pd$lineCol)), cex=1) }
# }
# dev.off()

## volcano plot
pdf("astrocytes_day77_volcano.pdf",h=6,w=6)
par(mar=c(5,6,5,2),cex.axis=1.5,cex.lab=1.4,cex.main=1.5)
plot(ffGene$logFC, -log10(ffGene$P.Value),
	 xlim=c(-15,15),
	 xlab="log2(FoldChange)", ylab="-log10(p-value)",
	 pch=21, col="black",	 
	 bg=ifelse(ffGene$bonf.P.Val<0.05, "orangered", NA) ) 
# abline(v=c(seq(-15,15, by=5)), lty=2, col="grey")

plot(ffGene$logFC, -log10(ffGene$P.Value),
	 xlim=c(-15,15),
	 xlab="log2(FoldChange)", ylab="-log10(p-value)",
	 pch=21, col="black",	 
	 bg=ifelse(ffGene$adj.P.Val<0.05, "orangered", NA) )
legend("topleft", "FDR < 0.05", pch=16, col="orangered", cex=.8)
# abline(v=c(seq(-15,15, by=5)), lty=2, col="grey")
dev.off()



################
# associations #
################

# gene set associations
library(clusterProfiler)
library(org.Hs.eg.db)

identical(ffGene$gencodeID, geneMap$gencodeID) ##check order

# split genes into up/down regulated, dropping grey
moduleGeneList = geneMap$EntrezID[ffGene$adj.P.Val < 0.05]
moduleGeneList = moduleGeneList[!is.na(moduleGeneList)]

# split genes into condition lists
up = geneMap$EntrezID[ffGene$adj.P.Val < 0.05 & ffGene$FC=="UP"]
down = geneMap$EntrezID[ffGene$adj.P.Val < 0.05 & ffGene$FC=="DOWN"]
moduleGeneList = list(up, down)
moduleGeneList = lapply(moduleGeneList, function(x) x[!is.na(x)])
names(moduleGeneList) = c("Up-regulated","Down-regulated")

## set universe of expressed genes
geneUniverse = as.character(geneMap$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

############################## 
## run enrichment analysis ###
##############################

## adjusted
goBP <- compareCluster(moduleGeneList, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Hs.eg.db,
                ont = "BP", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = .001,
				readable= TRUE)
goMF <- compareCluster(moduleGeneList, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Hs.eg.db,
                ont = "MF", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = .001,
				readable= TRUE)
goCC <- compareCluster(moduleGeneList, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Hs.eg.db,
                ont = "CC", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = .001,
				readable= TRUE)
kegg <- compareCluster(moduleGeneList, fun = "enrichKEGG",
                universe = geneUniverse,  pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = .001)
save(goBP,goMF,goCC,kegg, file="compareClusters_astro.rda")

load("compareClusters_astro.rda")

# ### combine
# geneSetStats_Adj = rbind(as.data.frame(goBP_Adj), as.data.frame(goMF_Adj),
	# as.data.frame(goCC_Adj), as.data.frame(kegg_Adj))
# geneSetStats_Adj$Ontology = rep(c("BP", "MF","CC", "KEGG"), 
	# times = c(nrow(goBP_Adj), nrow(goMF_Adj), nrow(goCC_Adj), nrow(kegg_Adj)))
# geneSetStats_Adj$Cluster = paste0("Adj_", geneSetStats_Adj$Cluster)

pdf("astro_up_down_enrichments_stemcell.pdf",h=5,w=9)
dotplot(goBP, includeAll="TRUE")
dotplot(goMF, includeAll="TRUE")
dotplot(goCC, includeAll="TRUE")
dotplot(kegg, includeAll="TRUE")
dev.off()




