###
library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(genefilter)
library(RColorBrewer)

MAINDIR = "/dcl01/lieber/ajaffe/lab/libd_stem_timecourse"


########################
### read in pheno
library(readxl)
## read in biological phenotype info
pheno = read_excel("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/prep_samples/MASTER_AZ_RNA-seq_10NOV_2016.xlsx")
## rename columns, some have changed
names(pheno)[c(9:10,12,16:28,31,43)] = c("LibConstruct",
	"DateToAJ", "JC_Comment", "NumReads", "SeqNotes",
	"NumMapped", "mapRate", "NumProperMap", "properRate",
	"NumUnmapped", "unmapRate",
	"riboMapped", "riboRate", "mitoMapped", "mitoRate",
	"alignNote", "LINE","Fluidigm")
## drop old tophat alignment info
pheno = pheno[,c(4:15, 29:37, 39:ncol(pheno))]
## make RNA Number and flowcell the unique ID
pheno$SampleID = paste0(pheno$RNA_NO, "_", pheno$Flowcell)


########################
### full samples
load("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/rse_gene_AZpilot_jan6.hg38_n506.Rdata")
load("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/rse_exon_AZpilot_jan6.hg38_n506.Rdata")
load("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/rse_jx_AZpilot_jan6.hg38_n506.Rdata")
pd = colData(rse_gene)

pd$SampleID = paste0(ss(pd$SAMPLE_ID, "_", 1), "_", 
	ss(pd$SAMPLE_ID, "_", 2))
	
########################
### combine info
pd = merge(pheno, pd, all=FALSE, by="SampleID")


########################
### organize
pd = pd[-which(pd$DAY == "NA" | pd$LINE =="165-B-8X"),]

pd$DAY = as.factor(paste0("day_",pd$DAY))
pd$DAY = factor(pd$DAY, levels(pd$DAY)[c(2,4,7,10,1,3,5,6,8,9)])

pd$Class[pd$Class=="ORF_MANIPULATION"] = "OE"
pd$Class[pd$Class=="shRNA_MANIPULATION"] = "SH"
pd$Class[pd$Class %in% c("Internal Control","Naked genomes")] = "TC"

pd$LINE[pd$SPECIES=="RAT"] = "RAT"

pd$COND = pd$CONDITION
pd$COND[pd$COND=="DORSAL(1)"] = "DORSAL"
pd$COND[pd$COND=="ACC_DORSAL(2)"] = "ACC_DORSAL"
pd$COND[grep("PLUS",pd$COND)] = "NEURONS+ASTROS"
pd$COND[pd$COND=="ASTROS_ALONE"] = "RAT_ASTROS"
pd$COND[pd$COND=="NEURONS_ALONE"] = "NEURONS"

pd$COND = as.factor(pd$COND)
pd$COND = factor(pd$COND,levels(pd$COND)[c(7,2,1,5,8,3,4,6)])

pd$LABEL = paste0(pd$COND,"_", pd$Class)
pd$LABEL_EXP = as.factor(pd$LABEL)
pd$LABEL_DAY = as.factor(pd$LABEL)

## group by OE/SH/timecourse
pd$LABEL_EXP = factor(pd$LABEL_EXP, levels(pd$LABEL_EXP)[c(4,1,12,5,2,8,13,11,3,9,14,6,7,10)])
## group by condition
pd$LABEL_DAY = factor(pd$LABEL_DAY, levels(pd$LABEL_DAY)[c(11,4:5,1:3,8:9,12:14,6:7,10)])



########################
### deconvolution


### more specific target - PTEN / NRXN1

pd = pd[!pd$COND %in% c("NEURONS","NEURONS+ASTROS","RAT_ASTROS"), ]
pd = pd[pd$Manipulation %in% c("NRXN1","PTEN","NO","shRNA_CNT"),]
pd$Manipulation[pd$Manipulation == "NO"] = "TC"

# Colors by Donor
pal = c(colorRampPalette(brewer.pal(10,"Paired")[1:2])(4), 
		colorRampPalette(brewer.pal(10,"Paired")[3:4])(3), 
		brewer.pal(10,"Paired")[9:10], brewer.pal(10,"Paired")[5:8])
pd$lineCol = pal[as.numeric(factor(pd$LINE))]


pd$LABEL_SH = paste0(pd$Manipulation,"_",pd$COND)
pd$LABEL_SH = as.factor(pd$LABEL_SH)
pd$LABEL_SH = factor(pd$LABEL_SH, levels(pd$LABEL_SH)[c(2,1,3,4,6,5,7,8,10,9,11,12,15,13,14,16)])

pd$LABEL_SH2 = factor(pd$LABEL_SH, levels(pd$LABEL_SH)[c(2,1,3,4,6,5,7,8,10,9,11,12,15,13,14,16)])

## time order
pd$LABEL_SH3 = factor(pd$LABEL_SH, levels(pd$LABEL_SH)[c(13,9,1,5,14,10,2,6,15,11,3,7,16,12,4,8)])

## filter
rse_gene = rse_gene[,match(pd$SAMPLE_ID, rse_gene$SAMPLE_ID)]
rse_exon = rse_exon[,match(pd$SAMPLE_ID, rse_exon$SAMPLE_ID)]
rse_jx = rse_jx[,match(pd$SAMPLE_ID, rse_jx$SAMPLE_ID)]
rowData(rse_jx)$Length = 100




####### PTEN plots

ptenInd = which(pd$Manipulation!="NRXN1" & !pd$EXPERIMENT %in% c("SH1","SH2","SH4","SH8") & 
				!pd$Donor%in%c("66","3"))
pd_pten = pd[ptenInd,]
# dropInd = grep("shRNA", pd_pten$LABEL_SH)
# pd_pten = pd_pten[-dropInd,]
pd_pten$LABEL_SH3 = droplevels(pd_pten$LABEL_SH3)
pd_pten$LABEL_SH3 = factor(pd_pten$LABEL_SH3, levels(pd_pten$LABEL_SH3)[c(3,2,6,5,1,4,9,8,7,12,11,10)])


rse_gene_pten = rse_gene[,ptenInd]
rse_exon_pten = rse_exon[,ptenInd]
rse_jx_pten = rse_jx[which(as.data.frame(rowRanges(rse_jx))$seqnames=="chr10"),ptenInd]
# rse_gene_pten = rse_gene_pten[,-dropInd]
# rse_exon_pten = rse_exon_pten[,-dropInd]
# rse_jx_pten = rse_jx_pten[,-dropInd]

yExprsG = log2(recount::getRPKM(rse_gene_pten,length_var="Length")+1)
yExprsE = log2(recount::getRPKM(rse_exon_pten,length_var="Length")+1)
	rse = rse_jx_pten
	length_var = "Length"
	mapped = colSums(as.data.frame(assays(rse)$counts))
	bg <- matrix(mapped, ncol = ncol(rse), nrow = nrow(rse),
		byrow = TRUE)
	len <- rowData(rse)[, length_var]
	wid <- matrix(len, nrow = nrow(rse), ncol = ncol(rse), byrow = FALSE)
	jRpkm = as.data.frame(assays(rse)$counts)/(wid/1000)/(bg/1e+06)
yExprsJ = log2(jRpkm+1)


# j1 = "chr10:89,653,867-89,690,802(+)"
# 		87,894,110 - 87,931,045

#j2 = "chr10:89,717,777-89,720,650(+)"
#		87,958,020 - 87,960,893

pdf("pten_features.pdf", h=5,w=6)
par(mar=c(11,4,2,2))

geneInd = which(rowData(rse_gene_pten)$Symbol=="PTEN")
boxplot(yExprsG[geneInd,] ~ pd_pten$LABEL_SH3,
        cex.axis=.8,cex.lab=1,las=3, main = "PTEN",
        ylab="log2(Exprs+1)", outline=FALSE )
points(yExprsG[geneInd,] ~ jitter(as.numeric(factor(pd_pten$LABEL_SH3))),
       pch=as.numeric(factor(pd_pten$SPECIES))+20,
       bg=pd_pten$lineCol, cex=1)
abline(v=c(2.5,6.5,9.5), col="grey", lty=2)   

exonInd = which(rowRanges(rse_exon_pten)$Symbol == "PTEN")
for (i in 1:length(exonInd)) {
	exoni = exonInd[i]
	boxplot(yExprsE[exoni,] ~ pd_pten$LABEL_SH3,
			cex.axis=.8,cex.lab=1,las=3, main=rownames(rowData(rse_exon_pten))[exoni],
			ylab="log2(Exprs+1)", outline=FALSE )
	points(yExprsE[exoni,] ~ jitter(as.numeric(factor(pd_pten$LABEL_SH3))),
		   pch=as.numeric(factor(pd_pten$SPECIES))+20,
		   bg=pd_pten$lineCol, cex=1)
	abline(v=c(2.5,6.5,9.5), col="grey", lty=2)   
}

jInd = which(rowRanges(rse_jx_pten)$newGeneSymbol == "PTEN")
for (i in 1:length(jInd)) {
	jxi = jInd[i]
	boxplot(as.numeric(yExprsJ[jxi,]) ~ pd_pten$LABEL_SH3,
			cex.axis=.8,cex.lab=1,las=3, main=rownames(rowData(rse_jx_pten))[jxi],
			ylab="log2(Exprs+1)", outline=FALSE )
	points(as.numeric(yExprsJ[jxi,]) ~ jitter(as.numeric(factor(pd_pten$LABEL_SH3))),
		   pch=as.numeric(factor(pd_pten$SPECIES))+20,
		   bg=pd_pten$lineCol, cex=1)
	abline(v=c(2.5,6.5,9.5), col="grey", lty=2)   
}

dev.off()



####### NRXN1 plots

nrxnInd = which(pd$Manipulation != "PTEN" & !pd$EXPERIMENT %in% c("SH3") & pd$Donor!="3")
pd_nrxn = pd[nrxnInd,]
# dropInd = grep("shRNA", pd_nrxn$LABEL_SH3)
# pd_nrxn = pd_nrxn[-dropInd,]
pd_nrxn$LABEL_SH3 = droplevels(pd_nrxn$LABEL_SH3)
pd_nrxn$LABEL_SH3 = factor(pd_nrxn$LABEL_SH3, levels(pd_nrxn$LABEL_SH3)[c(3,2,6,5,1,4,9,8,7,12,11,10)])

rse_gene_nrxn = rse_gene[,nrxnInd]
rse_exon_nrxn = rse_exon[,nrxnInd]
rse_jx_nrxn = rse_jx[which(as.data.frame(rowRanges(rse_jx))$seqnames=="chr2"),nrxnInd]
# rse_gene_nrxn = rse_gene_nrxn[,-dropInd]
# rse_exon_nrxn = rse_exon_nrxn[,-dropInd]
# rse_jx_nrxn = rse_jx_nrxn[,-dropInd]

yExprsG = log2(recount::getRPKM(rse_gene_nrxn,length_var="Length")+1)
yExprsE = log2(recount::getRPKM(rse_exon_nrxn,length_var="Length")+1)
	rse = rse_jx_nrxn
	length_var = "Length"
	mapped = colSums(as.data.frame(assays(rse)$counts))
	bg <- matrix(mapped, ncol = ncol(rse), nrow = nrow(rse),
		byrow = TRUE)
	len <- rowData(rse)[, length_var]
	wid <- matrix(len, nrow = nrow(rse), ncol = ncol(rse), byrow = FALSE)
	jRpkm = as.data.frame(assays(rse)$counts)/(wid/1000)/(bg/1e+06)
yExprsJ = log2(jRpkm+1)



## first exon

pdf("nrxn_features.pdf", h=5,w=6)
par(mar=c(11,4,2,2))
palette(brewer.pal(8,"Dark2"))

geneInd = which(rowData(rse_gene_nrxn)$Symbol=="NRXN1")
boxplot(yExprsG[geneInd,] ~ pd_nrxn$LABEL_SH3,
        cex.axis=.8,cex.lab=.8,las=3, main = "NRXN1",
        ylab="log2(Exprs+1)", outline=FALSE )
points(yExprsG[geneInd,] ~ jitter(as.numeric(factor(pd_nrxn$LABEL_SH3))),
       pch=as.numeric(factor(pd_nrxn$SPECIES))+20,
       bg=pd_nrxn$lineCol, cex=1)
abline(v=c(2.5,6.5,9.5), col="grey", lty=2)   

exonInd = which(rowRanges(rse_exon_nrxn)$Symbol == "NRXN1")
for (i in 1:length(exonInd)) {
	exoni = exonInd[i]
	boxplot(yExprsE[exoni,] ~ pd_nrxn$LABEL_SH3, ylim=c(min(yExprsE[exoni,]),max(yExprsE[exoni,])),
			cex.axis=.8,cex.lab=.8,las=3, main=rownames(rowData(rse_exon_nrxn))[exoni],
			ylab="log2(Exprs+1)", outline=FALSE )
	points(yExprsE[exoni,] ~ jitter(as.numeric(factor(pd_nrxn$LABEL_SH3))),
		   pch=as.numeric(factor(pd_nrxn$SPECIES))+20,
		   bg=pd_nrxn$lineCol, cex=1)
	abline(v=c(2.5,6.5,9.5), col="grey", lty=2)   
}

jInd = which(rowRanges(rse_jx_nrxn)$newGeneSymbol == "NRXN1")
for (i in 1:length(jInd)) {
	jxi = jInd[i]
	boxplot(as.numeric(yExprsJ[jxi,]) ~ pd_nrxn$LABEL_SH3,
			cex.axis=.8,cex.lab=.8,las=3, main=rownames(rowData(rse_jx_nrxn))[jxi],
			ylab="log2(Exprs+1)", outline=FALSE )
	points(as.numeric(yExprsJ[jxi,]) ~ jitter(as.numeric(factor(pd_nrxn$LABEL_SH3))),
		   pch=as.numeric(factor(pd_nrxn$SPECIES))+20,
		   bg=pd_nrxn$lineCol, cex=1)
	abline(v=c(2.5,6.5,9.5), col="grey", lty=2)   
}

dev.off()


















########################
### deconvolution

###############################################################
############## cell type ######################################

pal = brewer.pal(10, "Set3")
pal[2] = "#8B795E" # darken NPCs

load("cell_type/singleCell_iPSC_quake_coefEsts_calibration_Zscale.rda")

yExprsTC_Z = scale(yExprsTC[rownames(coefEsts),])

# project
stemPropEsts = minfi:::projectCellType(yExprsTC_Z,coefEsts)
stemPropEstsScaled = prop.table(stemPropEsts,1)

## line plot
gIndexes_tc = splitit(pd$LABEL_SH)
stemPropEsts_groupMeans = sapply(gIndexes_tc, 
	function(ii) colMeans(stemPropEsts[ii,]))
stemPropEsts_groupSEs = sapply(gIndexes_tc, 
	function(ii) apply(stemPropEsts[ii,],2,function(x) sd(x)/sqrt(length(x))))

round(100*stemPropEsts_groupMeans,3)


pdf("cell_type/lineplots/libd_stemcell_timecourse_cellTypeDecon_n257.pdf",h=5,w=7)
palette(pal)
par(mar=c(10,6,2,2),cex.axis=1.2,cex.lab=1.6)
plot(stemPropEsts_groupMeans[1,], type="p",xaxt="n",
	lty=1, pch=19,lwd=2, ylim = c(0,1),col=1,xlab="",
	ylab = "Class Proportion",cex=1.3)
axis(1,at=seq(along=gIndexes_tc),levels(pd$LABEL_SH), las=3,cex.axis=.78)
segments(x0=seq(along=gIndexes_tc), x1=seq(along=gIndexes_tc), col=1,lwd=2, 
	y0=stemPropEsts_groupMeans[1,] - 2*stemPropEsts_groupSEs[1,], 
	y1=stemPropEsts_groupMeans[1,] + 2*stemPropEsts_groupSEs[1,])
for(i in 1:10) {
	lines(stemPropEsts_groupMeans[i,1:4], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=seq(along=gIndexes_tc[1:4]), x1=seq(along=gIndexes_tc[1:4]), col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,1:4] - 2*stemPropEsts_groupSEs[i,1:4], 
	y1=stemPropEsts_groupMeans[i,1:4] + 2*stemPropEsts_groupSEs[i,1:4])
}
for(i in 1:10) {
	lines(x=5:8, y=stemPropEsts_groupMeans[i,5:8], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=5:8, x1=5:8, col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,5:8] - 2*stemPropEsts_groupSEs[i,5:8], 
	y1=stemPropEsts_groupMeans[i,5:8] + 2*stemPropEsts_groupSEs[i,5:8])
}
for(i in 1:10) {
	lines(x=9:12, y=stemPropEsts_groupMeans[i,9:12], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=9:12, x1=9:12, col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,9:12] - 2*stemPropEsts_groupSEs[i,9:12], 
	y1=stemPropEsts_groupMeans[i,9:12] + 2*stemPropEsts_groupSEs[i,9:12])
}
for(i in 1:10) {
	lines(x=13:16, y=stemPropEsts_groupMeans[i,13:16], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=13:16, x1=13:16, col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,13:16] - 2*stemPropEsts_groupSEs[i,13:16], 
	y1=stemPropEsts_groupMeans[i,13:16] + 2*stemPropEsts_groupSEs[i,13:16])
}
dev.off()




### by time
## line plot
gIndexes_tc = splitit(pd$LABEL_SH3)
stemPropEsts_groupMeans = sapply(gIndexes_tc, 
	function(ii) colMeans(stemPropEsts[ii,]))
stemPropEsts_groupSEs = sapply(gIndexes_tc, 
	function(ii) apply(stemPropEsts[ii,],2,function(x) sd(x)/sqrt(length(x))))

round(100*stemPropEsts_groupMeans,3)


pdf("cell_type/lineplots/libd_stemcell_timecourse_cellTypeDecon_n257_bytime.pdf",h=5,w=7)
palette(pal)
par(mar=c(10,6,2,2),cex.axis=1.2,cex.lab=1.6)
plot(stemPropEsts_groupMeans[1,], type="p",xaxt="n",
	lty=1, pch=19,lwd=2, ylim = c(0,1),col=1,xlab="",
	ylab = "Class Proportion",cex=1.3)
axis(1,at=seq(along=gIndexes_tc),levels(pd$LABEL_SH3), las=3,cex.axis=.78)
segments(x0=seq(along=gIndexes_tc), x1=seq(along=gIndexes_tc), col=1,lwd=2, 
	y0=stemPropEsts_groupMeans[1,] - 2*stemPropEsts_groupSEs[1,], 
	y1=stemPropEsts_groupMeans[1,] + 2*stemPropEsts_groupSEs[1,])
for(i in 1:10) {
	lines(x=1, y=stemPropEsts_groupMeans[i,1], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=1, x1=1, col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,1] - 2*stemPropEsts_groupSEs[i,1], 
	y1=stemPropEsts_groupMeans[i,1] + 2*stemPropEsts_groupSEs[i,1])
}
for(i in 1:10) {
	lines(x=2:4, y=stemPropEsts_groupMeans[i,2:4], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=2:4, x1=2:4, col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,2:4] - 2*stemPropEsts_groupSEs[i,2:4], 
	y1=stemPropEsts_groupMeans[i,2:4] + 2*stemPropEsts_groupSEs[i,2:4])
}
for(i in 1:10) {
	lines(x=5:8, y=stemPropEsts_groupMeans[i,5:8], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=5:8, x1=5:8, col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,5:8] - 2*stemPropEsts_groupSEs[i,5:8], 
	y1=stemPropEsts_groupMeans[i,5:8] + 2*stemPropEsts_groupSEs[i,5:8])
}
for(i in 1:10) {
	lines(x=9:12, y=stemPropEsts_groupMeans[i,9:12], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=9:12, x1=9:12, col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,9:12] - 2*stemPropEsts_groupSEs[i,9:12], 
	y1=stemPropEsts_groupMeans[i,9:12] + 2*stemPropEsts_groupSEs[i,9:12])
}
for(i in 1:10) {
	lines(x=13:16, y=stemPropEsts_groupMeans[i,13:16], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=13:16, x1=13:16, col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,13:16] - 2*stemPropEsts_groupSEs[i,13:16], 
	y1=stemPropEsts_groupMeans[i,13:16] + 2*stemPropEsts_groupSEs[i,13:16])
}
dev.off()





###############################################################
############## brain stage ####################################

pal = brewer.pal(8,"Dark2")

load("brain_stage/iPSC_BrainSpan_coefEsts_calibration_Zscore.rda")

yExprsTC_Z = scale(yExprsTC[rownames(coefEsts),])

# project
stemPropEsts = minfi:::projectCellType(yExprsTC_Z,coefEsts)
stemPropEstsScaled = prop.table(stemPropEsts,1)

## line plot
gIndexes_tc = splitit(pd$LABEL_SH)
stemPropEsts_groupMeans = sapply(gIndexes_tc, 
	function(ii) colMeans(stemPropEsts[ii,]))
stemPropEsts_groupSEs = sapply(gIndexes_tc, 
	function(ii) apply(stemPropEsts[ii,],2,function(x) sd(x)/sqrt(length(x))))

round(100*stemPropEsts_groupMeans,3)


pdf("brain_stage/lineplots/libd_stemcell_timecourse_brainStageDecon_n257.pdf",h=5,w=7)
palette(pal)
par(mar=c(10,6,2,2),cex.axis=1.2,cex.lab=1.6)
plot(stemPropEsts_groupMeans[1,], type="p",xaxt="n",
	lty=1, pch=19,lwd=2, ylim = c(0,1),col=1,xlab="",
	ylab = "Class Proportion",cex=1.3)
axis(1,at=seq(along=gIndexes_tc),levels(pd$LABEL_SH), las=3,cex.axis=.78)
segments(x0=seq(along=gIndexes_tc), x1=seq(along=gIndexes_tc), col=1,lwd=2, 
	y0=stemPropEsts_groupMeans[1,] - 2*stemPropEsts_groupSEs[1,], 
	y1=stemPropEsts_groupMeans[1,] + 2*stemPropEsts_groupSEs[1,])
for(i in 1:8) {
	lines(stemPropEsts_groupMeans[i,1:4], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=seq(along=gIndexes_tc[1:4]), x1=seq(along=gIndexes_tc[1:4]), col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,1:4] - 2*stemPropEsts_groupSEs[i,1:4], 
	y1=stemPropEsts_groupMeans[i,1:4] + 2*stemPropEsts_groupSEs[i,1:4])
}
for(i in 1:8) {
	lines(x=5:8, y=stemPropEsts_groupMeans[i,5:8], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=5:8, x1=5:8, col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,5:8] - 2*stemPropEsts_groupSEs[i,5:8], 
	y1=stemPropEsts_groupMeans[i,5:8] + 2*stemPropEsts_groupSEs[i,5:8])
}
for(i in 1:8) {
	lines(x=9:12, y=stemPropEsts_groupMeans[i,9:12], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=9:12, x1=9:12, col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,9:12] - 2*stemPropEsts_groupSEs[i,9:12], 
	y1=stemPropEsts_groupMeans[i,9:12] + 2*stemPropEsts_groupSEs[i,9:12])
}
for(i in 1:8) {
	lines(x=13:16, y=stemPropEsts_groupMeans[i,13:16], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=13:16, x1=13:16, col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,13:16] - 2*stemPropEsts_groupSEs[i,13:16], 
	y1=stemPropEsts_groupMeans[i,13:16] + 2*stemPropEsts_groupSEs[i,13:16])
}
dev.off()


### by time
## line plot
gIndexes_tc = splitit(pd$LABEL_SH3)
stemPropEsts_groupMeans = sapply(gIndexes_tc, 
	function(ii) colMeans(stemPropEsts[ii,]))
stemPropEsts_groupSEs = sapply(gIndexes_tc, 
	function(ii) apply(stemPropEsts[ii,],2,function(x) sd(x)/sqrt(length(x))))

round(100*stemPropEsts_groupMeans,3)


pdf("brain_stage/lineplots/libd_stemcell_timecourse_brainStageDecon_n257_bytime.pdf",h=5,w=7)
palette(pal)
par(mar=c(10,6,2,2),cex.axis=1.2,cex.lab=1.6)
plot(stemPropEsts_groupMeans[1,], type="p",xaxt="n",
	lty=1, pch=19,lwd=2, ylim = c(0,1),col=1,xlab="",
	ylab = "Class Proportion",cex=1.3)
axis(1,at=seq(along=gIndexes_tc),levels(pd$LABEL_SH3), las=3,cex.axis=.78)
segments(x0=seq(along=gIndexes_tc), x1=seq(along=gIndexes_tc), col=1,lwd=2, 
	y0=stemPropEsts_groupMeans[1,] - 2*stemPropEsts_groupSEs[1,], 
	y1=stemPropEsts_groupMeans[1,] + 2*stemPropEsts_groupSEs[1,])
for(i in 1:8) {
	lines(x=1, y=stemPropEsts_groupMeans[i,1], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=1, x1=1, col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,1] - 2*stemPropEsts_groupSEs[i,1], 
	y1=stemPropEsts_groupMeans[i,1] + 2*stemPropEsts_groupSEs[i,1])
}
for(i in 1:8) {
	lines(x=2:4, y=stemPropEsts_groupMeans[i,2:4], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=2:4, x1=2:4, col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,2:4] - 2*stemPropEsts_groupSEs[i,2:4], 
	y1=stemPropEsts_groupMeans[i,2:4] + 2*stemPropEsts_groupSEs[i,2:4])
}
for(i in 1:8) {
	lines(x=5:8, y=stemPropEsts_groupMeans[i,5:8], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=5:8, x1=5:8, col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,5:8] - 2*stemPropEsts_groupSEs[i,5:8], 
	y1=stemPropEsts_groupMeans[i,5:8] + 2*stemPropEsts_groupSEs[i,5:8])
}
for(i in 1:8) {
	lines(x=9:12, y=stemPropEsts_groupMeans[i,9:12], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=9:12, x1=9:12, col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,9:12] - 2*stemPropEsts_groupSEs[i,9:12], 
	y1=stemPropEsts_groupMeans[i,9:12] + 2*stemPropEsts_groupSEs[i,9:12])
}
for(i in 1:8) {
	lines(x=13:16, y=stemPropEsts_groupMeans[i,13:16], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=13:16, x1=13:16, col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,13:16] - 2*stemPropEsts_groupSEs[i,13:16], 
	y1=stemPropEsts_groupMeans[i,13:16] + 2*stemPropEsts_groupSEs[i,13:16])
}
dev.off()







##### Run voom on TC vs knockdowns
library(limma)
library(edgeR)

geneMap = as.data.frame(rowRanges(rse_gene))
geneRpkm = recount::getRPKM(rse_gene, "Length")
geneMap$meanExprs = rowMeans(geneRpkm)
gIndex = which(geneMap$meanExprs > 0.1)
rse_gene = rse_gene[gIndex,]

geneMap = as.data.frame(rowRanges(rse_gene))
geneCounts = assays(rse_gene)$counts

pd$Manip = as.factor(pd$Manipulation)
pd$Manip = relevel(pd$Manip, ref="TC")


############ ACC_DORSAL
subInd = which(pd$COND == "ACC_DORSAL")
pdSub = pd[subInd,]
geneCountsSub = geneCounts[,subInd]
## using contrasts
mod = model.matrix(~1 + Manip + LINE + totalAssignedGene, data=pd)

dgeGene = DGEList(counts = geneCounts, genes = geneMap)
dgeGene = calcNormFactors(dgeGene)
vGene = voom(dgeGene, mod, plot=FALSE)
fitGene = lmFit(vGene)
ebGene = ebayes(fitGene)
ffGene = topTable(eBayes(fitGene), coef=2:4, n = nrow(dgeGene)) #by condition
ffGene$bonf.P.Val = p.adjust(ffGene$P.Value, "bonf")
ffGene = ffGene[rownames(ebGene$t),]

voom_ad = cbind(ffGene, ebGene$t[,2:4], ebGene$p.value[,2:4])
names(voom_ad)[24:26] = paste0("t_", names(voom_ad)[24:26])
names(voom_ad)[27:29] = paste0("pvalue_", names(voom_ad)[27:29])
voom_ad$fdr_ManipNRXN1 = p.adjust(voom_ad$pvalue_ManipNRXN1, "fdr")
voom_ad$fdr_ManipPTEN = p.adjust(voom_ad$pvalue_ManipPTEN, "fdr")
voom_ad$fdr_ManipshRNA_CNT = p.adjust(voom_ad$pvalue_ManipshRNA_CNT, "fdr")




############ NPC
subInd = which(pd$COND == "NPC")
pdSub = pd[subInd,]
geneCountsSub = geneCounts[,subInd]
## using contrasts
mod = model.matrix(~1 + Manip + LINE + totalAssignedGene, data=pd)

dgeGene = DGEList(counts = geneCounts, genes = geneMap)
dgeGene = calcNormFactors(dgeGene)
vGene = voom(dgeGene, mod, plot=FALSE)
fitGene = lmFit(vGene)
ebGene = ebayes(fitGene)
ffGene = topTable(eBayes(fitGene), coef=2:4, n = nrow(dgeGene)) #by condition
ffGene$bonf.P.Val = p.adjust(ffGene$P.Value, "bonf")
ffGene = ffGene[rownames(ebGene$t),]

voom_npc = cbind(ffGene, ebGene$t[,2:4], ebGene$p.value[,2:4])
names(voom_npc)[24:26] = paste0("t_", names(voom_npc)[24:26])
names(voom_npc)[27:29] = paste0("pvalue_", names(voom_npc)[27:29])
voom_npc$fdr_ManipNRXN1 = p.adjust(voom_npc$pvalue_ManipNRXN1, "fdr")
voom_npc$fdr_ManipPTEN = p.adjust(voom_npc$pvalue_ManipPTEN, "fdr")
voom_npc$fdr_ManipshRNA_CNT = p.adjust(voom_npc$pvalue_ManipshRNA_CNT, "fdr")



############ ROSETTE
subInd = which(pd$COND == "ROSETTE")
pdSub = pd[subInd,]
geneCountsSub = geneCounts[,subInd]
## using contrasts
mod = model.matrix(~1 + Manip + LINE + totalAssignedGene, data=pd)

dgeGene = DGEList(counts = geneCounts, genes = geneMap)
dgeGene = calcNormFactors(dgeGene)
vGene = voom(dgeGene, mod, plot=FALSE)
fitGene = lmFit(vGene)
ebGene = ebayes(fitGene)
ffGene = topTable(eBayes(fitGene), coef=2:4, n = nrow(dgeGene)) #by condition
ffGene$bonf.P.Val = p.adjust(ffGene$P.Value, "bonf")
ffGene = ffGene[rownames(ebGene$t),]

voom_rose = cbind(ffGene, ebGene$t[,2:4], ebGene$p.value[,2:4])
names(voom_rose)[24:26] = paste0("t_", names(voom_rose)[24:26])
names(voom_rose)[27:29] = paste0("pvalue_", names(voom_rose)[27:29])
voom_rose$fdr_ManipNRXN1 = p.adjust(voom_rose$pvalue_ManipNRXN1, "fdr")
voom_rose$fdr_ManipPTEN = p.adjust(voom_rose$pvalue_ManipPTEN, "fdr")
voom_rose$fdr_ManipshRNA_CNT = p.adjust(voom_rose$pvalue_ManipshRNA_CNT, "fdr")







#### Cell type numbers
# > round(100*stemPropEsts_groupMeans,3)
                  # NRXN1_DORSAL NRXN1_ACC_DORSAL NRXN1_NPC NRXN1_ROSETTE PTEN_DORSAL PTEN_ACC_DORSAL PTEN_NPC PTEN_ROSETTE TC_RENEW TC_ACC_DORSAL TC_NPC TC_ROSETTE
# iPSC                    59.323           46.506    24.460        19.053      59.690          45.698   25.576       21.800   97.610        59.333 25.356     19.776
# NPC                     40.245           50.825    73.182        72.590      40.876          52.216   73.162       72.605    0.789        40.348 72.440     70.964
# Fetal_replicating        0.000            0.000     2.085         8.189       0.000           0.000    0.620        5.447    0.000         0.117  1.176      3.025
# Fetal_quiescent          1.384            2.317     2.244         7.537       0.315           2.008    0.713        8.332    2.220         1.531  1.405      6.441
# OPC                      0.000            0.000     0.000         0.000       0.000           0.000    0.000        0.000    0.000         0.000  0.000      0.000
# Neurons                  5.055            3.681     0.124         0.528       5.016           2.608    0.000        0.000    0.188         1.758  0.036      2.263
# Astrocytes               0.050            0.000     0.001         0.000       0.322           0.000    0.000        0.000    0.006         0.153  0.462      0.112
# Oligodendrocytes         0.000            0.000     0.000         0.000       0.000           0.000    0.000        0.000    0.038         0.002  0.000      0.000
# Microglia                0.000            0.000     0.000         0.000       0.000           0.000    0.000        0.000    0.000         0.000  0.000      0.000
# Endothelial              8.678           10.702    13.992         9.403       8.956          12.094   16.300       10.358    2.892         6.260 15.029     16.986


#### Brain stage numbers
# > round(100*stemPropEsts_groupMeans,3)
               # NRXN1_DORSAL NRXN1_ACC_DORSAL NRXN1_NPC NRXN1_ROSETTE PTEN_DORSAL PTEN_ACC_DORSAL PTEN_NPC PTEN_ROSETTE TC_RENEW TC_ACC_DORSAL TC_NPC TC_ROSETTE
# iPSC                 75.934           66.308    55.912        46.462      76.109          67.086   57.413       49.286   96.139        75.517 56.802     49.652
# NCX_EarlyFetal       24.509           31.260    36.193        32.087      24.460          29.003   32.881       33.128    0.457        20.259 33.815     31.999
# NCX_MidFetal          0.000            0.000     0.822        18.784       0.000           0.000    0.585       15.425    0.000         0.097  0.000      7.876
# NCX_LateFetal         0.000            0.602    10.015         4.227       0.000           2.707   11.720        5.239    0.000         1.196  9.848     12.531
# NCX_Infant            0.195            1.546     4.395         3.778       0.295           1.283    4.326        5.304    0.000         1.315  1.803      0.316
# NCX_Child             0.000            0.000     0.000         0.000       0.000           0.000    0.000        0.000    0.000         0.000  0.000      0.000
# NCX_Teens             0.000            0.113     0.000         0.000       0.000           0.000    0.000        0.000    0.000         0.076  0.000      0.000
# NCX_Adult             0.000            0.045     0.000         0.000       0.000           0.000    0.000        0.000    0.000         0.010  0.000      0.000







      0.000    0.000         0.000  0.000      0.000
# NCX_Teens             0.000            0.113     0.000         0.000       0.000           0.000    0.000        0.000    0.000         0.076  0.000      0.000
# NCX_Adult             0.000            0.045     0.000         0.000       0.000           0.000    0.000        0.000    0.000         0.010  0.000      0.000







