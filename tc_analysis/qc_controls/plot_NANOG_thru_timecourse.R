## load packages
library(jaffelab)
library(GenomicRanges)
library(SummarizedExperiment)
library(RColorBrewer)

MAINDIR = '/dcl01/lieber/ajaffe/lab/libd_stem_timecourse'
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n157.rda"))

## recount getRPKM version
getRPKM <- function(rse, length_var = 'bp_length', mapped_var = NULL) {
    mapped <- if(!is.null(mapped_var)) colData(rse)[, mapped_var] else colSums(assays(rse)$counts)      
    bg <- matrix(mapped, ncol = ncol(rse), nrow = nrow(rse), byrow = TRUE)   
    len <- if(!is.null(length_var)) rowData(rse)[, length_var] else width(rowRanges(rse))   
    wid <- matrix(len, nrow = nrow(rse), ncol = ncol(rse), byrow = FALSE)
    assays(rse)$counts / (wid/1000) / (bg/1e6)
}
geneRpkm = getRPKM(rse_gene, length_var="Length")

gMap = rowData(rse_gene)
pd = colData(rse_gene)
pd$REP = as.numeric(pd$BIO_REP)
pd$REP[is.na(pd$REP)] = 1
pd$DAY = ordered(pd$DAY, levels=sort(as.numeric(levels(as.factor(pd$DAY)))) )
pd$DIV = as.numeric(as.character(pd$DAY))

## Don't plot RENEW samples
keepInd = which(pd$CONDITION != "RENEW")
geneRpkm = geneRpkm[,keepInd]
pd = pd[keepInd,]

yExprs = as.matrix(log2(geneRpkm+1))
lineIndexes = splitit(paste0(pd$LINE,":",pd$REP))

# Colors by Donor
pal = c(colorRampPalette(brewer.pal(10,"Paired")[1:2])(5), 
		colorRampPalette(brewer.pal(10,"Paired")[3:4])(3), 
		brewer.pal(10,"Paired")[9:10], brewer.pal(10,"Paired")[5:8])
pd$lineCol = pal[as.numeric(factor(pd$LINE))]


## separate color legend
pdf("lineColoring.pdf")
plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
legend("topleft", levels(factor(pd$LINE)),
	pch = 15, col = pal, cex=1.8,bty="n",pt.cex=2.5)
dev.off()



##################################################
###### plots with only tc n=106 samples #######
library(jaffelab)
library(GenomicRanges)
library(rtracklayer)
library(SummarizedExperiment)
library(RColorBrewer)

MAINDIR = '/dcl01/lieber/ajaffe/lab/libd_stem_timecourse'
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n146.rda"))
pd = colData(rse_gene)
pd$DAY = as.factor(pd$DAY)
pd$DAY = factor(pd$DAY, levels(as.factor(pd$DAY))[c(2,4,6,9,1,3,5,7,8)])
gMap = rowData(rse_gene)

## recount getRPKM version
getRPKM <- function(rse, length_var = 'bp_length', mapped_var = NULL) {
    mapped <- if(!is.null(mapped_var)) colData(rse)[, mapped_var] else colSums(assays(rse)$counts)      
    bg <- matrix(mapped, ncol = ncol(rse), nrow = nrow(rse), byrow = TRUE)   
    len <- if(!is.null(length_var)) rowData(rse)[, length_var] else width(rowRanges(rse))   
    wid <- matrix(len, nrow = nrow(rse), ncol = ncol(rse), byrow = FALSE)
    assays(rse)$counts / (wid/1000) / (bg/1e6)
}
geneRpkm = getRPKM(rse_gene, length_var="Length")


## don't use technical replicates
bioInd = which(pd$Class == "Naked genomes")
pd = pd[bioInd,]
geneRpkm = geneRpkm[,bioInd]   # n=128

## don't use RENEW controls or NEURONS_ALONE in analyses
dropInd = which(pd$CONDITION %in% c("RENEW","NEURONS_ALONE"))
pd = pd[-dropInd,]
geneRpkm = geneRpkm[,-dropInd]   # n=106

yExprs = as.matrix(log2(geneRpkm+1))


# Colors by Donor
pal = c(colorRampPalette(brewer.pal(10,"Paired")[1:2])(4), 
		colorRampPalette(brewer.pal(10,"Paired")[3:4])(3), 
		brewer.pal(10,"Paired")[9:10], brewer.pal(10,"Paired")[5:8])
pd$lineCol = pal[as.numeric(factor(pd$LINE))]


#### List from Josh
# 1) early downregulated- POU5F1, NANOG
# 2) transiently early upregulated- HES5, NR2E1, DLL1, LGR5, SOX1
# 3) early upregulated and maintained- FOXG1, ZIC1
# 4) mid upregulated and maintained- NFIA, SLITRK3
#### List from http://docs.abcam.com/pdf/neuroscience/neural-lineage-markers-web.pdf
# DCX and TBR1 (immature neuron)
# RBFOX3 and MAP2 (mature neuron)
# SLC6A1 (GABAergic neuron)
# SLC18A3 (Cholinergic neuron)
# KCNJ6 (dopaminergic neuron)
# SLC17A6 (glutamatergic neurons)

genelist = c("POU5F1","NANOG","HES5","NR2E1","DLL1","LGR5","SOX1","FOXG1","ZIC1","NFIA","SLITRK3",
			"DCX","TBR1","RBFOX3","MAP2","SLC6A1","SLC18A3","KCNJ6","SLC17A6")	
label = c(rep("pluripotent",2), rep("neuronal precursor",5),
		rep("early upregulated and maintained",2), rep("mid upregulated and maintained",2),
		rep("immature neuron",2),rep("mature neuron",2),"GABAergic neuron","cholinergic neuron",
		"dopaminergic neuron","glutamatergic neuron")

pdf("CONTROLS_through_timecourse_n106_consistentAxis.pdf",h=7,w=7, useDingbats=FALSE)
par(mar=c(6,6,4,2),cex.axis=1.8,cex.lab=2,cex.main=2)
for (i in 1:length(genelist)) {
ind = which(gMap$Symbol==genelist[i])
plot(yExprs[ind,] ~ jitter(as.numeric(pd$DAY),.5), 
	pch = 21, bg=pd$lineCol,
	cex=1.75, xlab="Day",
	ylab="log2(Exprs + 1)", ylim=c(0,9),
	main = paste0(genelist[i], "\n", label[i]), xaxt="n")
axis(1, at=1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
if (i==1) {legend("topright", levels(factor(pd$LINE)),
	pch = 15, col = pal, cex=1.5,nc=2) }
}
dev.off()

pdf("CONTROLS_through_timecourse_n106_onepage.pdf",h=8,w=18, useDingbats=FALSE)
par(mfrow=c(2,4),mar=c(6,6,6,2),cex.axis=1.8,cex.lab=2,cex.main=2)
for (i in 1:length(genelist)) {
ind = which(gMap$Symbol==genelist[i])
plot(yExprs[ind,] ~ jitter(as.numeric(pd$DAY),.5), 
	pch = 21, bg=pd$lineCol,
	cex=2.5, xlab="Day",
	ylab="log2(Exprs + 1)", ylim=c(0,9),
	main = paste0(genelist[i], "\n", label[i]), xaxt="n")
axis(1, at=1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
if (i==2) {legend("topright", levels(factor(pd$LINE)),
	pch = 15, col = pal, cex=1.5,nc=2) }
}
dev.off()


### Paper figure
genelist = c("POU5F1","HES5","SLC17A6")
label = c("POU5F1/OCT3","HES5","SLC17A6/VGLUT2")
		
pdf("figure1_POU5F1_HES5_SLC17A6_n106.pdf",h=6,w=7, useDingbats=FALSE)
par(mar=c(6,6,4,2),cex.axis=1.8,cex.lab=2,cex.main=2)
for (i in 1:length(genelist)) {
ind = which(gMap$Symbol==genelist[i])
plot(yExprs[ind,] ~ as.numeric(pd$DAY), 
	pch = 21, bg=pd$lineCol,
	cex=2, xlab="Day",
	ylab="log2(Exprs + 1)", ylim=c(0,9),
	main = label[i], xaxt="n")
axis(1, at=1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
if (i==1) {legend("topright", levels(factor(pd$LINE)),
	pch = 15, col = pal, cex=1.3,nc=2) }
}
dev.off()



############################
## pca ########

pca1 = prcomp(t(log2(geneRpkm+1)))
pcaVars1 = getPcaVars(pca1)

pd$dayLabel = ifelse(pd$SPECIES=="HUMAN", paste0("Day: ", pd$DAY), pd$SPECIES)
pd$dayLabel = ordered(as.factor(pd$dayLabel), levels=levels(as.factor(pd$dayLabel))[c(2,4,5,6,1,3,7)])
levels(pd$dayLabel)[8] = "Rat Astros"

pd$dayLabel = paste0("Day: ", pd$DAY)
pd$dayLabel = ordered(as.factor(pd$dayLabel), levels=levels(as.factor(pd$dayLabel))[c(2,4,6,9,1,3,5,7,8)])

# Colors by Donor
daypal = c(brewer.pal(10,"Spectral")[c(1:4,6,7:10)])
pd$dayCol = daypal[as.numeric(factor(pd$dayLabel))]

pd$COND = as.factor(pd$CONDITION)
pd$COND = factor(pd$COND,levels(pd$COND)[c(1,3,4,2)])
levels(pd$COND) = c("Accelerated Dorsal","NPC","Rosette","Neurons on Rat Astros")

## PC 1 vs 2: Explains days / cell conditions
pdf("pca_log2Rpkm_PC1_2_n106.pdf", h=6,w=6)
par(mar=c(5,6,4,2),cex.axis=2,cex.lab=2,cex.main=2)
plot(pca1$x, pch = as.numeric(pd$COND)+20, bg=pd$dayCol,cex=2, main="Gene PCs",
     xlab=paste0("PC1: ", pcaVars1[1], "% Var Expl"),
     ylab=paste0("PC2: ", pcaVars1[2], "% Var Expl"))
legend("top", paste0(levels(factor(pd$dayLabel))),
       pch = 15, col = daypal, cex=.9)
legend("topright", paste0(levels(factor(pd$COND))),
       pch = c(16,15,18,17), cex=.85)

dev.off()





