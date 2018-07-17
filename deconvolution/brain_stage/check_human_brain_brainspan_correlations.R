#####

library(SummarizedExperiment)
library(recount)
library(jaffelab)
library(RColorBrewer)
library(minfi)
library(lattice)

MAINDIR = "/dcl01/lieber/ajaffe/lab/libd_stem_timecourse"

# ## human data from span
# load("/dcl01/lieber/ajaffe/lab/brainspan_analysis/rawCounts_full_brainspan_n607.rda", verbose=TRUE)
# load("/dcl01/lieber/ajaffe/lab/brainspan_analysis/full.brainspan_final.phenotype.rda",verbose=TRUE)
# pdSpan = cbind(fin_pdSpan, metrics[fin_pdSpan$lab,])
# pdSpan$wig = NULL
# geneCounts = geneCounts[,pdSpan$lab]
# rownames(pdSpan) = pdSpan$lab
# rm(exonCounts, exonMap, jCounts, jMap)

# ## Create gene,exon RangedSummarizedExperiment objects
# gr_genes <- GRanges(seqnames = geneMap$Chr,
    # IRanges(geneMap$Start, geneMap$End), strand = geneMap$Strand)
# names(gr_genes) <- rownames(geneMap)
# mcols(gr_genes) <- DataFrame(geneMap[, - which(colnames(geneMap) %in%
    # c('Chr', 'Start', 'End', 'Strand'))])
# rse_geneSpan <- SummarizedExperiment(assays = list('counts' = geneCounts),
    # rowRanges = gr_genes, colData = pdSpan)
	
# ### using AgeCat
# rse_geneSpan$AgeCat = "Adult"
# rse_geneSpan$AgeCat[rse_geneSpan$Age<20] = "Teens"
# rse_geneSpan$AgeCat[rse_geneSpan$Age<10] = "Child"
# rse_geneSpan$AgeCat[rse_geneSpan$Age<1] = "Infant"
# rse_geneSpan$AgeCat[rse_geneSpan$Age<0] = "Late Fetal"
# rse_geneSpan$AgeCat[rse_geneSpan$Age<(-.25)] = "Mid Fetal"
# rse_geneSpan$AgeCat[rse_geneSpan$Age<(-.5)] = "Early Fetal"
# rse_geneSpan$AgeCat = as.factor(rse_geneSpan$AgeCat)
# rse_geneSpan$AgeCat = factor(rse_geneSpan$AgeCat,levels(rse_geneSpan$AgeCat)[c(3,6,5,4,2,7,1)])

# ## brainspan NCX
# rse_geneSpan$RegionGroup = rse_geneSpan$Regioncode 
# rse_geneSpan$RegionGroup[rse_geneSpan$Regioncode %in% 
				# c("A1C", "DFC", "IPC","ITC", "M1C", "MFC",
				# "OFC", "S1C", "STC", "V1C", "VFC")] = "NCX"

# save(rse_geneSpan, file = "brainspan_gene_rse.rda")
load("brainspan_gene_rse.rda")

geneRpkmSpan = getRPKM(rse_geneSpan,length_var="Length")
yExprsSpan= log2(geneRpkmSpan+1)

################
## bring in scorecard data
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


############
## filter and merge

## just iPSC
rse_geneSC = rse_geneSC[,rse_geneSC$CONDITION == "iPSC"] #
yExprsSC = log2(getRPKM(rse_geneSC, "Length")+1)


rse_geneSpan$AgeRegion = paste0(rse_geneSpan$RegionGroup,"_", 
					gsub(" ", "", rse_geneSpan$AgeCat))
rse_geneSpan$AgeRegion = factor(rse_geneSpan$AgeRegion,
	paste0(rep(unique(rse_geneSpan$RegionGroup),each= 7), "_", 
		rep(gsub(" ", "", levels(rse_geneSpan$AgeCat)), times = 7)))
gIndexes = splitit(rse_geneSpan$AgeRegion)
gIndexes = gIndexes[lengths(gIndexes) > 1]
meanExprsSpan = sapply(gIndexes, function(ii) rowMeans(yExprsSpan[,ii]))
meanExprsSem = rowMeans(yExprsSC)
meanExprs = cbind(iPSC = meanExprsSem, meanExprsSpan)

proteinIndex = which(rowData(rse_geneSC)$gene_type == "protein_coding")

##################################
### project into brainseq ########
##################################

theSeq = seq(-1,1,by=0.05)
my.col <- colorRampPalette(c("blue","white","red"))(length(theSeq))

### validate using BrainSeq
load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseGene_merged_n732.rda")
rse_geneDLPFC = rse_gene[,rse_gene$Dx == "Control" & 
	rse_gene$Race %in% c("AA","CAUC")]
rse_geneDLPFC$ageGroup = cut(rse_geneDLPFC$Age, c(-1,-0.45, 0,1,10,20,100))
levels(rse_geneDLPFC$ageGroup) = c("EarlyFetal","MidFetal",
	"Infant","Child","Teens","Adult")

rse_geneDLPFC = rse_geneDLPFC[,order(rse_geneDLPFC$Age)]
yExprsDLPFC = log2(getRPKM(rse_geneDLPFC,length_var="Length")+1)

meanExprsDLPFC = sapply(splitit(rse_geneDLPFC$ageGroup), function(ii) rowMeans(yExprsDLPFC[,ii]))

ccDlpfc = cor(meanExprsDLPFC, meanExprs)
ccDlpfcProt = cor(meanExprsDLPFC[proteinIndex,], meanExprs[proteinIndex,])
signif(ccDlpfc,3)
signif(ccDlpfcProt,3)
print(levelplot(ccDlpfc, aspect = "fill", at = theSeq,pretty=TRUE,
	panel = panel.levelplot.raster, col.regions = my.col,
		scales=list(x=list(rot=90)),	ylab = "Cell Type", xlab = ""))
print(levelplot(ccDlpfcProt, aspect = "fill", at = theSeq,pretty=TRUE,
	panel = panel.levelplot.raster, col.regions = my.col,
		scales=list(x=list(rot=90)),	ylab = "Cell Type", xlab = ""))

###################
## based on stem ##
###################
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n146.rda"))
## don't use technical replicates
bioInd = which(colData(rse_gene)$Class == "Naked genomes")
rse_gene = rse_gene[,bioInd]  # n=128
colData(rse_gene) = colData(rse_gene)
colData(rse_gene)$DAY = paste0("day_",colData(rse_gene)$DAY)
colData(rse_gene)$experiment = "lieber"
rse_geneTC = rse_gene

## clean CONDITION
rse_geneTC$COND = rse_geneTC$CONDITION
rse_geneTC$COND[grep("ACC_DORSAL",rse_geneTC$COND)] = "ACC_DORSAL"
rse_geneTC$COND[grep("PLUS",rse_geneTC$COND)] = "NEURONS+ASTROS"
rse_geneTC$COND = as.factor(rse_geneTC$COND)
rse_geneTC$COND = factor(rse_geneTC$COND,levels(rse_geneTC$COND)[c(5,1,4,6,2,3)])

## rpkm function
getRPKM = recount::getRPKM

geneRpkmTC = getRPKM(rse_geneTC, "Length")
yRpkmTC = log2(geneRpkmTC+1)

meanExprsStem = sapply(splitit(rse_geneTC$COND), function(ii) rowMeans(yRpkmTC[,ii]))
ccStem = cor(meanExprsStem, meanExprs)
print(levelplot(ccStem, aspect = "fill", at = theSeq,pretty=TRUE,
	panel = panel.levelplot.raster, col.regions = my.col,
		scales=list(x=list(rot=90)),	ylab = "Cell Type", xlab = ""))

pdf("stemCell_correlationToBrainspan.pdf",w=10)
ccStemProtein = cor(meanExprsStem[proteinIndex,], meanExprs[proteinIndex,])
print(levelplot(ccStemProtein, aspect = "fill", at = theSeq,pretty=TRUE,
	panel = panel.levelplot.raster, col.regions = my.col,
		ylab = "BrainSpan Sample", xlab = "Stem Cell Stage"))
dev.off()
