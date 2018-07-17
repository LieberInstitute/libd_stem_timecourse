####
library(SummarizedExperiment)
library(jaffelab)

pd = read.csv("phenoData.csv")

load("pairedend_n180/rse_gene_yeo_singlecell_pe_n180.Rdata")
rse_paired = rse_gene
load("singleend_n34/rse_gene_yeo_singlecell_se_n34.Rdata")
rse_single = rse_gene

## match up colData column names
colData(rse_single) = colData(rse_single)[,c(1:13,25:45)]
colData(rse_paired) = colData(rse_paired)[,c(1:13,36:58)]
colData(rse_single)$trimmed = NA
colData(rse_single)$concordMapRate = NA
colData(rse_single) = colData(rse_single)[,colnames(colData(rse_paired))]

## match up mcols (remove meanExprs)
mcols(rse_single) = mcols(rse_single)[,-8]
mcols(rse_paired) = mcols(rse_paired)[,-8]
## merge
rse_gene = cbind(rse_single, rse_paired)

## add phenoData.csv
rownames(pd) = pd$Run
pd = pd[colnames(rse_gene),]
colData(rse_gene) = cbind(pd, colData(rse_gene))

save(rse_gene, file="rse_gene_yeo_n214.Rdata")
