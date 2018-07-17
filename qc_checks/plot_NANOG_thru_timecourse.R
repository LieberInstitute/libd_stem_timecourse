## load packages
library(jaffelab)
library(GenomicRanges)
library(SummarizedExperiment)
library(RColorBrewer)

MAINDIR = '/dcl01/lieber/ajaffe/lab/libd_stem_timecourse'
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n157.rda"))

geneRpkm = getRPKM(rse_gene)
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

# ## separate color legend
# pdf("NANOGlineColoring.pdf")
# plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
# legend("topleft", levels(factor(pd$LINE)),
	# pch = 15, col = pal, cex=1.8,bty="n",pt.cex=2.5)
# dev.off()


#### Expression of NANOG gene across timecourse
genelist = "NANOG"	
NANOGind = which(gMap$Symbol %in% c("NANOG"))

## connect points from 165-B-8X
lineInd = which(pd$LINE=="165-B-8X")

## Day as factor
pdf("NANOG_through_timecourse.pdf",h=6,w=8)
par(mar=c(6,6,4,2),cex.axis=2,cex.lab=2,cex.main=2)
plot(as.numeric(pd$DAY), yExprs[NANOGind,], 
	pch = 21, bg=pd$lineCol,
	cex=2, xlab="Day",
	ylab="log2(Exprs + 1)",
	main = "NANOG", xaxt="n")
axis(1, at=1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
lines(as.numeric(pd$DAY)[lineInd], yExprs[NANOGind,lineInd], lwd=1.5, col=pal[5])
legend("topright", levels(factor(pd$LINE)),
	pch = 15, col = pal, cex=1)
dev.off()



#### List from Josh
# 1) early downregulated- POU5F1, NANOG
# 2) transiently early upregulated- HES5, NR2E1, DLL1, LGR5, SOX1
# 3) early upregulated and maintained- FOXG1, ZIC1
# 4) mid upregulated and maintained- NFIA, SLITRK3

genelist = c("POU5F1","NANOG","HES5","NR2E1","DLL1","LGR5","SOX1","FOXG1","ZIC1","NFIA","SLITRK3")	
label = c(rep("pluripotent",2), rep("neuronal precursor",5),
		rep("early upregulated and maintained",2), rep("mid upregulated and maintained",2))

pdf("CONTROLS_through_timecourse.pdf",h=6,w=8)
par(mar=c(6,6,4,2),cex.axis=2,cex.lab=2,cex.main=2)
for (i in 1:length(genelist)) {
ind = which(gMap$Symbol==genelist[i])
plot(as.numeric(pd$DAY), yExprs[ind,], 
	pch = 21, bg=pd$lineCol,
	cex=2, xlab="Day",
	ylab="log2(Exprs + 1)",
	main = paste0(genelist[i], "\n", label[i]), xaxt="n")
axis(1, at=1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
if (i==1) {legend("topright", levels(factor(pd$LINE)),
	pch = 15, col = pal, cex=1) }
}
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
geneRpkm = getRPKM(rse_gene)

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

pdf("CONTROLS_through_timecourse_n106_consistentAxis_aj.pdf",h=7,w=7, useDingbats=FALSE)
par(mar=c(6,6,4,2),cex.axis=1.8,cex.lab=2,cex.main=2)
for (i in 1:length(genelist)) {
ind = which(gMap$Symbol==genelist[i])
plot(yExprs[ind,] ~ jitter(as.numeric(pd$DAY),.5), 
	pch = 21, bg=pd$lineCol,
	cex=1.75, xlab="Day",
	ylab="log2(Exprs + 1)", ylim=c(0,7.5),
	main = paste0(genelist[i], "\n", label[i]), xaxt="n")
axis(1, at=1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
if (i==1) {legend("topright", levels(factor(pd$LINE)),
	pch = 15, col = pal, cex=1.5,nc=2) }
}
dev.off()









