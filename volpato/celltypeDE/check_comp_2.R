##
library(jaffelab)
library(minfi)
library(SummarizedExperiment)
library(recount)
library(limma)
library(edgeR)
library(RColorBrewer)

########################
### read in pheno
pheno = read.csv("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/volpato/volpato_pheno_bulk.csv", 
				stringsAsFactors=FALSE)
rownames(pheno) = pheno$SAMPLE_ID


########################
### full samples
load("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/volpato/bulk/rse_gene_volpato_bulk_n57.Rdata")
pd = colData(rse_gene)
pd = cbind(pheno[pd$SAMPLE_ID,], pd)

#get RPKM
yExprsTC = log2(recount::getRPKM(rse_gene,length_var="Length")+1)


pd$cell_line = as.factor(pd$cell_line)
pd$site = as.factor(pd$site)
pd$time_point = as.factor(pd$time_point)


#######################
##   deconvolution   ##
#######################

## cell type 
load("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/deconvolution/cell_type/singleCell_iPSC_quake_coefEsts_calibration_Zscale.rda")

yExprsTC_Z = scale(yExprsTC[rownames(coefEsts),])

# project
cellPropEsts = minfi:::projectCellType(yExprsTC_Z,coefEsts)
cellPropEsts = as.data.frame(cellPropEsts)

## merge
pd = cbind(pd, cellPropEsts)




#####################
## do DE modeling ###
#####################

## filter low expression
geneRpkm = recount::getRPKM(rse_gene,length_var="Length")
gIndex = rowMeans(geneRpkm) > 0.2
rse_gene = rse_gene[gIndex,]

geneCounts = assays(rse_gene)$counts
geneMap = as.data.frame(rowRanges(rse_gene))

#### voom ####
dgeGene = DGEList(counts = geneCounts, genes = geneMap)
dgeList = calcNormFactors(dgeGene)


####### model 1 - only site
mod0 = model.matrix(~ site, data=pd)
voomList = voom(dgeList, mod0, plot = FALSE)
geneFit = lmFit(voomList)
statList1_lab0 = topTable(eBayes(geneFit), coef=2:5, n = nrow(dgeGene))[,-15]
statList1_lab0$bonf.P.Val = p.adjust(statList1_lab0$P.Val, "bonferroni")
sum(statList1_lab0$bonf.P.Val < 0.05)
# [1] 7616
 
####### model 2 - only neuron
mod0 = model.matrix(~ Neurons, data=pd)
voomList = voom(dgeList, mod0, plot = FALSE)
geneFit = lmFit(voomList)
statList2_neur0 = topTable(eBayes(geneFit), coef=2, n = nrow(dgeGene))[,-15]
statList2_neur0$bonf.P.Val = p.adjust(statList2_neur0$P.Val, "bonferroni")
sum(statList2_neur0$bonf.P.Val < 0.05)
# [1] 2814

 
####### model 3 - site + neuron
## adjustment
mod = model.matrix(~ site + Neurons, data=pd)
voomList = voom(dgeList, mod, plot = FALSE)
geneFit = lmFit(voomList)
# Lab stats
statList3_lab = topTable(eBayes(geneFit), coef=2:5, n = nrow(dgeGene))[,-15]
statList3_lab$bonf.P.Val = p.adjust(statList3_lab$P.Val, "bonferroni")
sum(statList3_lab$bonf.P.Val < 0.05)
# [1] 5316

# Neuron stats
statList3_neur = topTable(eBayes(geneFit), coef=6, n = nrow(dgeGene))[,-15]
statList3_neur$bonf.P.Val = p.adjust(statList3_neur$P.Val, "bonferroni")
sum(statList3_neur$bonf.P.Val < 0.05)
# [1] 594



######################
## interaction model #
######################


####### model 4 - site + neuron + interaction

modInt = model.matrix(~ site*Neurons, data=pd)

voomList = voom(dgeList, modInt, plot = FALSE)
geneFit = lmFit(voomList)

# Lab stats
statList4_lab = topTable(eBayes(geneFit), coef=2:5, n = nrow(dgeGene))[,-15]
statList4_lab$bonf.P.Val = p.adjust(statList4_lab$P.Val, "bonferroni")
sum(statList4_lab$bonf.P.Val < 0.05)
# [1] 1835

# Neuron stats
statList4_neur = topTable(eBayes(geneFit), coef=6, n = nrow(dgeGene))[,-15]
statList4_neur$bonf.P.Val = p.adjust(statList4_neur$P.Val, "bonferroni")
sum(statList4_neur$bonf.P.Val < 0.05)
# [1] 52

# Interaction stats
statList4_int = topTable(eBayes(geneFit), coef=7:10, n = nrow(dgeGene))[,-15]
statList4_int$bonf.P.Val = p.adjust(statList4_int$P.Val, "bonferroni")
sum(statList4_int$bonf.P.Val < 0.05)
# [1] 1099





### model 1
sum(statList1_lab0$bonf.P.Val < 0.05)
# [1] 7616

### model 2
sum(statList2_neur0$bonf.P.Val < 0.05)
# [1] 2814

### model 3
sum(statList3_lab$bonf.P.Val < 0.05)
# [1] 5316
sum(statList3_neur$bonf.P.Val < 0.05)
# [1] 594

### model 4
sum(statList4_lab$bonf.P.Val < 0.05)
# [1] 1835
sum(statList4_neur$bonf.P.Val < 0.05)
# [1] 52
sum(statList4_int$bonf.P.Val < 0.05)
# [1] 1099



table(statList1_lab0$bonf.P.Val < 0.05, statList4_lab$bonf.P.Val < 0.05)







statListMain = topTable(eBayes(geneFit), coef=2, n = nrow(dgeGene))[,-15]
statListInt = topTable(eBayes(geneFit), coef=8, n = nrow(dgeGene))[,-15]
table(statListMain$adj.P.Val < 0.01)
# FALSE  TRUE
# 19413    22
table(statListInt$adj.P.Val < 0.01)
# FALSE
# 19435














# are these 22 in the previous significant lists
sigEns = c(statList0$gencodeID[statList0$adj.P.Val<0.01], statList$gencodeID[statList$adj.P.Val<0.01])
table(statListMain$gencodeID[statListMain$adj.P.Val<0.01] %in% sigEns)	
# 17 of 22


## merge
colnames(statListMain) = paste0(colnames(statListMain), "_main")
colnames(statListInt) = paste0(colnames(statListInt), "_int")
szStatsEx = cbind(statListMain, statListInt)

sigStats = szStatsEx[rowSums(szStatsEx[,grep("adj", colnames(szStatsEx))] < 0.01) > 0,]

sigStats = sigStats[order(sigStats$P.Value_main),!grepl("Ave|B", colnames(sigStats))]
write.csv(sigStats, file = "volpato_cellModeling_sigStats.csv")

#######################
## show some plots?? ##

geneExprs = as.matrix(log2(recount::getRPKM(rse_geneLab, "Length")+1))
pdLab$NeurQuad = cut(pdLab$Neurons, breaks = quantile(pdLab$Neurons), include=TRUE,labels=FALSE)
pdLab$DxQuad = paste0(pdLab$NeurQuad, ":", pdLab$site)

pdf("example_boxplots_of_interaction.pdf",width=10)
par(mar=c(5,6,4,2), cex.axis=1.8,cex.lab=1.8,cex.main=1.8)
palette(c("blue", "red"))
for(i in 1:nrow(sigStats)) {
	boxplot(geneExprs[rownames(sigStats)[i],] ~ pdLab$DxQuad,
		main = paste0(sigStats$Symbol[i]), ylab= "log2(RPKM+1)",
		xlab = "Neuron RNA Fraction Quantile & Site", # col=rep(1:2,times=4),
		outline=FALSE, ylim = range(geneExprs[rownames(sigStats)[i],]))
	points(geneExprs[rownames(sigStats)[i],] ~ 
		jitter(as.numeric(factor(pdLab$DxQuad)), amount=0.15), 
		pch = 21, bg=factor(ss(pdLab$DxQuad, ":", 2)))
	abline(v = c(2.5,3.5,5.5),lty=2,cex=1.8)
	legend("topleft", paste0(c("Main", "Int"),
		" p=", signif(sigStats[i,c(17,35)],3)),nc=1,cex=1.8)
}
dev.off()

