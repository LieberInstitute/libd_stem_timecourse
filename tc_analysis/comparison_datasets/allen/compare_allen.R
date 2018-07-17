###
library(minfi)
library(jaffelab)
library(SummarizedExperiment)
library(genefilter)
library(edgeR)
library(preprocessCore)
library(RColorBrewer)

######################
## load allen data
load("/dcl01/ajaffe/data/lab/singleCell/YaoAllen/rse_gene_allen_singleCell_hg38_n4556.Rdata")
pdAllen = colData(rse_gene)
geneCountsAllen = assays(rse_gene)$counts
geneMapAllen = as.data.frame(rowRanges(rse_gene))

## recount getRPKM version
getRPKM <- function(rse, length_var = 'bp_length', mapped_var = NULL) {
    mapped <- if(!is.null(mapped_var)) colData(rse)[, mapped_var] else colSums(assays(rse)$counts)      
    bg <- matrix(mapped, ncol = ncol(rse), nrow = nrow(rse), byrow = TRUE)   
    len <- if(!is.null(length_var)) rowData(rse)[, length_var] else width(rowRanges(rse))   
    wid <- matrix(len, nrow = nrow(rse), ncol = ncol(rse), byrow = FALSE)
    assays(rse)$counts / (wid/1000) / (bg/1e6)
}
geneRpkmAllen = getRPKM(rse_gene, length_var="Length")

pdAllen$proto = ss(as.character(pdAllen$title),"_",3)

pd1 = pdAllen[which(pdAllen$proto=="smart-seq"),]
pd2 = pdAllen[which(pdAllen$proto!="smart-seq" | is.na(pdAllen$proto)),]

# ## only use smart-seq protocol
# smartInd = grep("smart-seq",pdAllen$title)
# pdAllen = pdAllen[smartInd,]
# geneCountsAllen = geneCountsAllen[,smartInd] # n=1846
# geneRpkmAllen = geneRpkmAllen[,smartInd]
## remove low alignment rates
dropInd = which(pdAllen$overallMapRate<0.8 | pdAllen$totalAssignedGene<0.4 | pdAllen$mitoRate>0.1)
pdAllen = pdAllen[-dropInd,]
geneCountsAllen = geneCountsAllen[,-dropInd] # n=3807
geneRpkmAllen = geneRpkmAllen[,-dropInd]



#######################
## load timecourse data
MAINDIR = '/dcl01/lieber/ajaffe/lab/libd_stem_timecourse'
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n146.rda"))
geneCounts = assays(rse_gene)$counts
pd = colData(rse_gene)

## don't use technical replicates
bioInd = which(pd$Class == "Naked genomes")
pd = pd[bioInd,]
geneCounts = geneCounts[,bioInd]   # n=128
## don't use RENEW controls or NEURONS_ALONE in analyses
dropInd = which(pd$CONDITION %in% c("RENEW","NEURONS_ALONE"))
pd = pd[-dropInd,]
geneCounts = geneCounts[,-dropInd]   # n=106



yExprsAllen = log2(geneRpkmAllen+1)


#######################
#######  compare  #####

# maybe take a peek at the allen single cell data at some of our most time-course associated DE genes

load(file.path(MAINDIR,"tc_analysis/rda/voom_tc_genes.rda"))  # voomStatsC
voomStatsC = voomStatsC[order(voomStatsC$P.Value),]

## plot Allen for top 100 LIBD genes

pdf("allen_top_libd_tc_genes_8.pdf",h=6,w=6)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(8,"Spectral"))	
for (i in 1:10) {
genei = rownames(voomStatsC)[i]
 plot(jitter(as.numeric(as.factor(pdAllen$DIV))), yExprsAllen[genei,] , xaxt="n", 
	   pch = 21, bg=as.numeric(as.factor(pdAllen$DIV)),
       cex=2,xlab="Day",
       ylab="log2(Exprs + 1)",
       main = paste0(voomStatsC$Symbol[i], "\n", voomStatsC$gencodeID[i]) )
 axis(1, at=1:length(levels(as.factor(pdAllen$DIV))), labels = levels(as.factor(pdAllen$DIV)))

}
dev.off()


## which Allen genes are most different?
geneMapAllen$t_p = NA
for (i in 1:nrow(geneRpkmAllen)) {
	geneMapAllen$t_p[i] = t.test(geneRpkmAllen[i,] ~ as.factor(pdAllen$DIV))$p.value
	if (i%%1000 == 0) { print(i) }
}
save(geneMapAllen, file="../rda/allen_ttests.rda")

## which Allen genes are most different?
geneMapAllen$lm_p = NA
for (i in 1:nrow(geneRpkmAllen)) {
	geneMapAllen$lm_p[i] = lm(yExprsAllen[i,] ~ pdAllen$DIV)$coef[2]
	if (i%%1000 == 0) { print(i) }
}
save(geneMapAllen, file="../rda/allen_lmtests.rda")



geneMapAllen = geneMapAllen[order(geneMapAllen$lm_p),]

pdf("allen_tc_associated_d.pdf",h=6,w=6)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(8,"Spectral"))	
for (i in 1:5) {
genei = rownames(geneMapAllen)[i]
sm.density.compare(yExprsAllen[genei,], as.factor(pdAllen$DIV), lty=1)
legend("topright", levels(as.factor(pdAllen$DIV)), fill=2+(0:nlevels(as.factor(pdAllen$DIV))))
 # plot(jitter(as.numeric(as.factor(pdAllen$DIV))), yExprsAllen[genei,] , xaxt="n",
	   # pch = 21, bg=as.numeric(as.factor(pdAllen$DIV)),
       # cex=2,xlab="Day",
       # ylab="log2(Exprs + 1)",
       # main = paste0(geneMapAllen$Symbol[i], "\n", geneMapAllen$gencodeID[i]) )
 # axis(1, at=1:length(levels(as.factor(pdAllen$DIV))), labels = levels(as.factor(pdAllen$DIV)))

}
dev.off()




