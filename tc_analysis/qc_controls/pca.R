###

library(SummarizedExperiment)
library(recount)
library(jaffelab)
library(RColorBrewer)

## libd time course
MAINDIR = "/dcl01/lieber/ajaffe/lab/libd_stem_timecourse"
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n146.rda"))
## don't use technical replicates
bioInd = which(colData(rse_gene)$Class == "Naked genomes")
rse_gene = rse_gene[,bioInd]  # n=128
colData(rse_gene) = colData(rse_gene)[,c(24,26,23,15:17,20,6)]
colData(rse_gene)$DAY = paste0("day_",colData(rse_gene)$DAY)
colData(rse_gene)$experiment = "lieber"
rse_geneTC = rse_gene

## put conditions in order
rse_geneTC$CONDITION[grep("ACC_DORS",rse_geneTC$CONDITION)] = "ACC_DORSAL"
rse_geneTC$CONDITION = as.factor(rse_geneTC$CONDITION)
rse_geneTC$CONDITION = factor(rse_geneTC$CONDITION, levels = levels(rse_geneTC$CONDITION)[c(5,1,4,6,2,3)])

## put days in order
rse_geneTC$DIV = as.factor(rse_geneTC$DAY)
rse_geneTC$DIV = factor(rse_geneTC$DIV, levels = levels(rse_geneTC$DIV)[c(2,4,6,9,1,3,5,7,8)])

## rpkm function
getRPKM = recount::getRPKM
# ## recount getRPKM version
# getRPKM <- function(rse, length_var = 'bp_length', mapped_var = NULL) {
    # mapped <- if(!is.null(mapped_var)) colData(rse)[, mapped_var] else colSums(assays(rse)$counts)      
    # bg <- matrix(mapped, ncol = ncol(rse), nrow = nrow(rse), byrow = TRUE)   
    # len <- if(!is.null(length_var)) rowData(rse)[, length_var] else width(rowRanges(rse))   
    # wid <- matrix(len, nrow = nrow(rse), ncol = ncol(rse), byrow = FALSE)
    # assays(rse)$counts / (wid/1000) / (bg/1e6)
# }

#####################
## PCA based on LIBD ##
#####################

geneRpkmTC = getRPKM(rse_geneTC, "Length")
yRpkmTC = log2(geneRpkmTC+1)
pcaTC = prcomp(t(yRpkmTC))
pcaTC_Vars = getPcaVars(pcaTC)


## put conditions in order
rse_geneTC$COND = as.character(rse_geneTC$CONDITION)
rse_geneTC$COND = ifelse(rse_geneTC$COND == "ACC_DORSAL", paste0(rse_geneTC$COND,": ", rse_geneTC$DAY), rse_geneTC$COND)
rse_geneTC$COND = as.factor(rse_geneTC$COND)
rse_geneTC$COND = factor(rse_geneTC$COND, levels = levels(rse_geneTC$COND)[c(8,1:4,7,9,5,6)])



pdf("pca_log2Rpkm_n128.pdf", h=6,w=7)
par(mar=c(5,6,4,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(9,"Spectral"))

plot(pcaTC$x, bg = as.numeric(factor(rse_geneTC$CONDITION))+3, cex=2, main="Gene PCs",
	pch=ifelse(rse_geneTC$CONDITION=="ACC_DORSAL",as.numeric(factor(rse_geneTC$COND))+19,21),
	xlab=paste0("PC1: ", pcaTC_Vars[1], "% Var Expl"),
    ylab=paste0("PC2: ", pcaTC_Vars[2], "% Var Expl"))
legend("bottom", paste0(levels(factor(rse_geneTC$COND))), col=c(4,rep("black",4),6:9), pt.bg=5,
       pch = c(16,21:24,rep(16,4)), cex=.9)
	   
plot(pcaTC$x, bg = as.numeric(factor(rse_geneTC$CONDITION))+3, pch=21, cex=2, main="Gene PCs",
	xlab=paste0("PC1: ", pcaTC_Vars[1], "% Var Expl"),
    ylab=paste0("PC2: ", pcaTC_Vars[2], "% Var Expl"))
legend("bottom", paste0(levels(factor(rse_geneTC$CONDITION))), col=4:9,
       pch = 16, cex=.9)



dev.off()





##### line/day checker matrix
library(gplots)

## make data frame of 1s and 0s
dat = table(rse_geneTC$LINE,rse_geneTC$DIV)
dat = ifelse(dat>0,1,0)
colnames(dat) = gsub("_", " ",colnames(dat))


datBefore = dat
## add 165-B-8X and dropped genotypes
datBefore = rbind(datBefore, "165-B-8X" = c(1,0,1,1,1,1,0,0,0))
datBefore = rbind(datBefore, "RAT" = c(0,0,0,0,0,1,1,1,1))
datBefore["90-A-10",c("day 15","day 21")] = 1		
datBefore = datBefore[c(1:4,14,5:13,15),]


pdf("sample_heatmap.pdf", h=6,w=6)
par(mar=c(1,1,1,1),cex.axis=2,cex.lab=2,cex.main=1.5)
##Before dropping
heatmap.2(datBefore, Rowv=FALSE,Colv=FALSE, scale="none", key=FALSE, dendrogram="none", trace="none", 
		margins=c(5,8), lwid=c(1,5), lhei=c(1,5), srtCol=60, main="Samples, n=165",
		colsep=0:9, rowsep=0:15, col=c("white","gray40"), sepcolor="black", sepwidth=c(0.02,0.02))
##After dropping
heatmap.2(dat, Rowv=FALSE,Colv=FALSE, scale="none", key=FALSE, dendrogram="none", trace="none", 
		margins=c(5,8), lwid=c(1,5), lhei=c(1,5), srtCol=60, main="Samples, n=128",
		colsep=0:9, rowsep=0:13, col=c("white","gray40"), sepcolor="black", sepwidth=c(0.02,0.02))
dev.off()
		
		
		
		
		
		
		
 