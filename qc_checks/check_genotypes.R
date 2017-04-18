######
library(readxl)
library(jaffelab)
library(readr)
library(parallel)
library(ggplot2)
library(jaffelab)
library(GenomicRanges)
library(rtracklayer)
library(SummarizedExperiment)
library(RColorBrewer)

MAINDIR = '/dcl01/lieber/ajaffe/lab/libd_stem_timecourse'
DATADIR = '/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell'

load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n157.rda"))
pd = colData(rse_gene)
pd = pd[,1:24]

## read in genotype text files
fn <- file.path(DATADIR, 'Genotypes', paste0(pd$SAMPLE_ID, '_calledVariants.txt'))
names(fn) = pd$SampleID
stopifnot(all(file.exists(fn)))

genoDat = lapply(fn, function(x) { read.delim(x, header=TRUE, as.is=TRUE) } )
genoDat = lapply(genoDat, function(x) { cbind(x, name=paste0(x$seq,"_",x$pos))  } )
snvDat = lapply(genoDat, function(x) { x[,c(5:8,11)] } )


dd = dist(snvDat)
hc = hclust(dd)

pdf("dendrogram_genotype.pdf",h=5,w=7)
par(mar=c(8,5,2,2))
palette(brewer.pal(8,"Dark2"))
myplclust(hc, lab.col=as.numeric(as.factor(pd$Donor)), xlab="",
          lab = paste0(pd$SampleID,"_",pd$Donor), main = "")
legend("topright", paste0(levels(factor(pd$Donor))),
       pch=15, col=1:8, cex=1)  				   
dev.off()




