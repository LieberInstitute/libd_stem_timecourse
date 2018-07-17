###
library(SummarizedExperiment)
library(recount)
library(minfi)
library(GEOquery)
library(edgeR)
library(limma)
library(jaffelab)

########################
# load tekin data ######
########################

load("preprocessed_human/rse_gene_tekin_human_n157.Rdata")
rse_gene$Species[rse_gene$Species == "Human and mouse"] = "Mix"

## get expression
yExprs = log2(getRPKM(rse_gene,length_var="Length")+1)

## figure out phenotypes
table(rse_gene$age_days, rse_gene$Species)

## figure out better fields
rse_gene$isSorted = grepl("sorted", rse_gene$strain_cell_type)

## geo phenotypes
geoData = getGEO("GSE111831")
geo_pheno = do.call("rbind", lapply(geoData, function(x) pData(x)[,c("geo_accession","title")]))

rse_gene$geo_title = as.character(geo_pheno$title[match(rse_gene$Sample_Name, 
					geo_pheno$geo_accession)])
rse_gene$geo_label = ss(rse_gene$geo_title, "_rep")

## rename
rse_gene_NGN2 = rse_gene

###############################
## check composition first ####
###############################

#######################
## load timecourse data
MAINDIR = '/dcl01/lieber/ajaffe/lab/libd_stem_timecourse'
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n146.rda"))
geneCounts = assays(rse_gene)$counts
pd = colData(rse_gene)
getRPKM = recount::getRPKM
geneRpkm = getRPKM(rse_gene, length_var="Length")

## don't use technical replicates
bioInd = which(pd$Class == "Naked genomes")
pd = pd[bioInd,]
geneRpkm = geneRpkm[,bioInd]
geneCounts = geneCounts[,bioInd]   # n=128


############################################################
## don't use RENEW controls or NEURONS_ALONE in analyses #### didn't drop these when calculating pca ####
dropInd = which(pd$CONDITION %in% c("RENEW","NEURONS_ALONE"))
pd = pd[-dropInd,]
geneRpkm = geneRpkm[,-dropInd]
geneCounts = geneCounts[,-dropInd]   # n=106

## clean CONDITION
pd$COND = pd$CONDITION
pd$COND = as.factor(pd$COND)
pd$COND = factor(pd$COND,levels(pd$COND)[c(1,3,4,2)])  ## put levels in order
pd$DAY = ordered(pd$DAY, levels=sort(as.numeric(levels(as.factor(pd$DAY)))) )

## drop low expression
geneIndex = which(rowMeans(geneRpkm) > 0.1)  #### didn't drop these when calculating pca ####
geneRpkm = geneRpkm[geneIndex,]
geneCounts = geneCounts[geneIndex,]
geneMap = rowData(rse_gene)[geneIndex,]  ## 25466 genes

# (don't drop RENEW in libd) # n=128
## clean CONDITION
pd$COND = pd$CONDITION
pd$COND[grep("ACC_DORSAL",pd$COND)] = "ACC_DORSAL"
pd$COND = as.factor(pd$COND)
pd$COND = factor(pd$COND,levels(pd$COND)[c(5,1,4,6,2,3)])  ## put levels in order
pd$DAY = ordered(pd$DAY, levels=sort(as.numeric(levels(as.factor(pd$DAY)))) )

######################
##  PCA ####
#################
geneRpkmNgn2 = getRPKM(rse_gene_NGN2, "Length")
geneRpkmNgn2 = geneRpkmNgn2[rownames(geneRpkm),]

pcaTC = prcomp(t(log2(geneRpkm+1)))
pcaTC_vars = getPcaVars(pcaTC)

geneExprsNgn2_Scaled = scale(t(log2(geneRpkmNgn2+1)), pcaTC$center, pcaTC$scale) 
genePCs_Ngn2_projected = geneExprsNgn2_Scaled %*% pcaTC$rotation 

## first PC
par(mar=c(5,30,1,1))
boxplot(genePCs_Ngn2_projected[,1] ~ rse_gene_NGN2$geo_label,horizontal=TRUE,las=1)


par(mar=c(5,6,4,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(9,"Spectral"))
plot(pcaTC$x, main="Gene PCs", cex=2, 
	pch = 21,
	bg = as.numeric(pd$COND)+3,
	xlab=paste0("PC1: ", pcaTC_vars[1], "% Var Expl"),
    ylab=paste0("PC2: ", pcaTC_vars[2], "% Var Expl"))
legend("top", paste0(levels(factor(pd$COND))), col=4:9,
       pch=16, cex=.9)
palette(brewer.pal(9,"Set1"))
points(genePCs_Ngn2_projected[,1:2], pch = 22, cex=2,
	bg = as.numeric(factor(rse_gene_NGN2$culture_condition)))
legend("topright", paste0("Day ",levels(factor(pdCort$Day))), col=1:9,
       pch=15, cex=.9)
