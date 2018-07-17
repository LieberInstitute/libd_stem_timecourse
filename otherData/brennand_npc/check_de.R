##
library(GenomicRanges)

## our data
load("/users/ajaffe/Lieber/Projects/RNAseq/SzControl_DE_paper/rdas/all_de_features.rda")
geneStats = outStats$Gene

## their data
neuron = read.delim("COS_DE_SZNeuron.tsv",as.is=TRUE,row.names=1)
npc = read.delim("COS_DE_SZNPC.tsv",as.is=TRUE,row.names=1)
npc = npc[rownames(neuron),]

mm = match(rownames(neuron), names(geneStats))
neuron = neuron[!is.na(mm),]
npc = npc[!is.na(mm),]
geneStats = geneStats[mm[!is.na(mm)]]
identical(rownames(neuron), rownames(npc))

## 
cor(neuron$t, as.data.frame(mcols(geneStats)[,grep("tstat", colnames(mcols(geneStats)))]))
cor(npc$t, as.data.frame(mcols(geneStats)[,grep("tstat", colnames(mcols(geneStats)))]))

chisq.test(table(neuron$P.Value < 0.05, geneStats$pval_adj < 0.05))
chisq.test(table(neuron$P.Value < 0.01, geneStats$pval_qsva < 0.01))

table(geneStats$fdr_qsva < 0.1)
mean(neuron$P.Value[geneStats$fdr_qsva < 0.1] < 0.05, na.rm=TRUE)