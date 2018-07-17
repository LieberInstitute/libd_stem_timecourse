###
library(minfi)
library(jaffelab)
library(SummarizedExperiment)
library(genefilter)
library(edgeR)
library(preprocessCore)
library(RColorBrewer)
library(pheatmap)

MAINDIR = "/dcl01/lieber/ajaffe/lab/libd_stem_timecourse"

######################
## load stemcell data
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n146.rda"))
## don't use technical replicates
bioInd = which(colData(rse_gene)$Class == "Naked genomes")
rse_gene = rse_gene[,bioInd]  # n=128
colData(rse_gene) = colData(rse_gene)[,c(24,26,23,15,17,20,6)]
colData(rse_gene)$DAY = paste0("day_",colData(rse_gene)$DAY)
colData(rse_gene)$experiment = "lieber"
rse_geneTC = rse_gene

######################
## load ScoreCard data
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

rowData(rse_geneTC)$Symbol = rowData(rse_geneSC)$Symbol
rowData(rse_geneTC) = rowData(rse_geneTC)[,1:5]
rowData(rse_geneSC) = rowData(rse_geneSC)[,1:5]

######################
## combine
rse_gene = cbind(rse_geneSC,rse_geneTC)

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
pd$cond = as.character(pd$CONDITION)
pd$cond[grep("ACC_DORSAL",pd$cond)] = "ACC_DORSAL"
pd$cond = as.factor(pd$cond)
pd$cond = factor(pd$cond, levels = levels(pd$cond)[c(2,4,3,8,1,7,9,5,6)])


######################
## PCA
pca1 = prcomp(t(log2(geneRpkm+1)))
pcaVars1 = getPcaVars(pca1)

## PC 1 vs 2: Explains cell conditions (1) / experiment (2)
pdf("scorecard_pca_combined_data.pdf", h=6,w=6)
par(mar=c(5,6,4,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(9,"Spectral"))
plot(pca1$x, main="Gene PCs", cex=2,
	 pch=as.numeric(as.factor(pd$experiment))+20, bg=as.numeric(as.factor(pd$cond)), 
     xlab=paste0("PC1: ", pcaVars1[1], "% Var Expl"),
     ylab=paste0("PC2: ", pcaVars1[2], "% Var Expl"))
legend("bottomright", paste0(levels(factor(pd$cond))), col=1:9,
       pch = c(15,15,15,rep(16,6)), cex=.7)

dev.off()



# ######################
# ## Heatmap

# dd = dist(t(log2(geneRpkm+1)))
# dd = as.matrix(dd)

# rownames(dd) = colnames(dd) = paste0(pd$experiment,"_",pd$cond,"_",pd$DAY)

# pdf("scorecard_pheatmap.pdf",h=30,w=30)
# pheatmap(dd, 
		# cluster_rows=T, 
		# cluster_cols=T,
		# color = colorRampPalette(rev(brewer.pal(n = 5, name ="Blues")))(20))
# dev.off()




#### subset of 96 genes
pd$label = paste0(pd$experiment,"_",pd$cond,"_",pd$DAY)
pd$label[pd$cond=="RENEW"] = paste0(pd$experiment[pd$cond=="RENEW"],"_",pd$cond[pd$cond=="RENEW"])

geneMap = as.data.frame(rowRanges(rse_gene))
genelist = scan("scorecard_genes.txt", what="characters", sep="\n")
keepInd = which(geneMap$Symbol %in% genelist)
subRpkm = geneRpkm[keepInd,]

ddSub = dist(t(log2(subRpkm+1)))
ddSub = as.matrix(ddSub)
rownames(ddSub) = colnames(ddSub) = pd$label


# pdf("scorecard_pheatmap_96genes.pdf",h=30,w=30)
# pheatmap(ddSub, 
		# cluster_rows=T, 
		# cluster_cols=T,
		# color = colorRampPalette(rev(brewer.pal(n = 5, name ="Blues")))(20))
# dev.off()

######################
## PCA
pca1 = prcomp(t(log2(subRpkm+1)))
pcaVars1 = getPcaVars(pca1)

## PC 1 vs 2: Explains cell conditions (1) / experiment (2)
pdf("scorecard_pca_96genes2.pdf", h=6,w=6)
par(mar=c(5,6,4,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(9,"Spectral"))
plot(pca1$x, main="Gene PCs", cex=2,
	 pch=as.numeric(as.factor(pd$experiment))+20, bg=as.numeric(as.factor(pd$cond)), 
     xlab=paste0("PC1: ", pcaVars1[1], "% Var Expl"),
     ylab=paste0("PC2: ", pcaVars1[2], "% Var Expl"))
legend("topright", paste0(levels(factor(pd$cond))), col=1:9,
       pch = c(15,15,15,rep(16,6)), cex=.8)

dev.off()



### plot just AD to look at individual days

pd$day = ss(as.character(pd$DAY),"_",2)
youngInd = which(pd$day %in% c(NA,2,4,6,9) & pd$cond!="Fibroblast")
pcaSub = pca1$x[youngInd,]
pdSub = pd[youngInd,]
pdSub$day = ss(as.character(pdSub$DAY),"_",2)
pdSub$day[pdSub$cond=="ACC_DORSAL"] = paste0("ACC_DORSAL_",as.character(pdSub$DAY[pdSub$cond=="ACC_DORSAL"]))
pdSub$day[pdSub$experiment=="scorecard"] = as.character(pdSub$cond)[pdSub$experiment=="scorecard"]
pdSub$day[pdSub$cond=="RENEW"] = "RENEW"
pdSub$day = as.factor(pdSub$day)
pdSub$day = factor(pdSub$day, levels = levels(pdSub$day)[c(5,6,7,1:4)])

## PC 1 vs 2: Explains cell conditions (1) / experiment (2)
pdf("scorecard_pca_96genes.pdf", h=6,w=6)
par(mar=c(5,6,4,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(7,"Spectral"))
plot(pcaSub, main="Gene PCs", cex=2, ylim=c(-10,13), xlim=c(-14,12),
	 pch=as.numeric(as.factor(pdSub$experiment))+20, bg=as.numeric(as.factor(pdSub$day)), 
     xlab=paste0("PC1: ", pcaVars1[1], "% Var Expl"),
     ylab=paste0("PC2: ", pcaVars1[2], "% Var Expl"))
legend("topright", paste0(levels(factor(pdSub$day))), col=1:9,
       pch = c(15,15,rep(16,4)), cex=1)

dev.off()




##############################################
## pca w/in scorecard and then project LIBD ##
##############################################

geneRpkmSC = getRPKM(rse_geneSC,length_var="Length")
pcaSC = prcomp(t(log2(geneRpkmSC+1)))
pcaSC_vars = getPcaVars(pcaSC)

geneRpkmTC = getRPKM(rse_geneTC,length_var="Length")
geneRpkmTC_Scaled = scale(t(log2(geneRpkmTC+1)), pcaSC$center, pcaSC$scale) 
genePCs_TC_projected = geneRpkmTC_Scaled %*% pcaSC$rotation 

### EMILY, CLEAN UP AND MAKE LIKE OTHER PLOTS

pdf("scorecard_pca_projected.pdf", h=6,w=7)
par(mar=c(5,6,4,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(9,"Spectral"))
plot(pcaSC$x, cex=2, main="Gene PCs", 
	pch = ifelse(pd$DX[pd$experiment=="scorecard"]=="isogenic",22,24),
	bg = pd$cond[pd$experiment == "scorecard"],
	xlab=paste0("PC1: ", pcaSC_vars[1], "% Var Expl"),
    ylab=paste0("PC2: ", pcaSC_vars[2], "% Var Expl"))
points(genePCs_TC_projected[,1:2], pch = 21, cex=2,
	bg = pd$cond[pd$experiment != "scorecard"])
legend("topleft", paste0(levels(factor(pd$cond))), col=1:9,
       pch = c(rep(15,3),rep(16,8)), cex=.86)
legend("top", c("Isogenic", "Nonisogenic"), col=1,
       pch = c(22,24), cex=.86)
dev.off()


#################################################
## regression calibration for identity scores ###
#################################################

### split up
ySC = as.matrix(log2(geneRpkmSC+1))
pdSC = pd[pd$experiment == "scorecard",]
pdSC$cond = factor(as.character(pdSC$cond))
y = as.matrix(log2(geneRpkm+1))

pluripotent = ifelse(pdSC$cond %in% c("iPSC", "ESC"), 1,0)
somatic = ifelse(! pdSC$cond %in% c("iPSC", "ESC"), 1,0)
mod = data.frame(model.matrix(~pluripotent+somatic - 1))

library(limma)
fit = lmFit(ySC, mod)
fit1 = lmFit(ySC, model.matrix(~pluripotent)) # for filter
eb1 = ebayes(fit1) # fit filter
ind = which(eb1$p[,2] < 1e-60) # only those different
Xmat = fit$coef[ind,]
Dmat = t(Xmat)%*%Xmat
guess = apply(y[ind,], 2, function(x)  solve(Dmat, t(Xmat) %*% x))[1,]
pd$Pluripotent = guess


pd$day = ss(as.character(pd$DAY),"_",2)
pd$cond2 = as.character(pd$cond)
pd$cond2[pd$cond2=="ACC_DORSAL"] = paste0("ACC_DORSAL day ",as.character(pd$day[pd$cond2=="ACC_DORSAL"]))
pd$cond2 = as.factor(pd$cond2)
pd$cond2 = factor(pd$cond2, levels = levels(pd$cond2)[c(5,7,6,11,1:4,10,12,8,9)])


## EMILY: make boxplot nice for paper, 
## 	overly jitter, color, etc
pdf("scorecard_boxplot_identityScores.pdf", h=6,w=8)
par(mar=c(12,6,2,2),cex.axis=1.2,cex.lab=2,cex.main=2)
palette(brewer.pal(9,"Spectral"))
boxplot(pd$Pluripotent ~ pd$cond2, las=2, outline=FALSE,
ylab="Pluripotent identity", xaxt="n")
points(pd$Pluripotent ~ jitter(as.numeric(pd$cond2)),
        pch=ifelse(pd$experiment=="scorecard",22,21), 
        bg=as.numeric(pd$cond), cex=2)
abline(v=3.5, col="gray", lty=2)
axis(1, at=1:12, labels = FALSE)
text(1:12, par("usr")[3]-.07,
     srt = 60, adj= 1, xpd = TRUE,
     labels = levels(pd$cond2), cex=1.2)
dev.off()

## numbers for text: 
tapply(pd$Pluripotent, pd$cond2, mean)
tapply(pd$Pluripotent, pd$cond2, sd)

summary(lm(Pluripotent ~ factor(cond2), data =pd, experiment != "scorecard"))
