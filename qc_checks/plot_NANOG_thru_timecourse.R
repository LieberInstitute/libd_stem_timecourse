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
label = c(rep("early downregulated",2), rep("transiently early upregulated",5),
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



