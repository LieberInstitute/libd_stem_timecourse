###

library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(genefilter)
library(RColorBrewer)

pal = brewer.pal(10, "Set3")
pal[2] = "#8B795E" # darken NPCs
 
MAINDIR = "/dcl01/lieber/ajaffe/lab/libd_stem_timecourse"

# load gene-level
load("/dcl01/ajaffe/data/lab/singleCell/Darmanis/rna-seq-pipeline_run2/rse_gene_Darmanis_scRNASeq_Darmanis_n466.Rdata")
rse_geneQuake = rse_gene
rm(rse_gene)

# phenotype
load("/dcl01/lieber/ajaffe/PublicData/Brain/Darmanis_Quake/rdas/phenotype_quake_processed.rda")
colData(rse_geneQuake) = cbind(colData(rse_geneQuake), 
	pd[match(colnames(rse_geneQuake), pd$RunName),2:14])

rse_geneQuake$Cell_type[rse_geneQuake$Cell_type == "hybrids"] = "hybrid"
## exclude hybrid
rse_geneQuake = rse_geneQuake[,rse_geneQuake$Cell_type != "hybrid"]

# get expression
yExprsQuake = log2(getRPKM(rse_geneQuake, "Length")+1)

## bring in YEO iPSC data
load(file.path(MAINDIR,"yeo_singleCell/rse_gene_yeo_n214.Rdata"))

rse_geneYeo = rse_gene[,rse_gene$cell_type %in% c("iPSC", "NPC") &
	rse_gene$sample_type=="Single_Cell"]
yExprsYeo = log2(getRPKM(rse_geneYeo, "Length")+1)

##########
## merge #

## merge w/ scorecard
yExprs_Merge = cbind(yExprsYeo, yExprsQuake)	
group = c(as.character(rse_geneYeo$cell_type), 
	as.character(rse_geneQuake$Cell_type))
	
group = factor(group,levels =c("iPSC", "NPC","Fetal_replicating",
	"Fetal_quiescent", "OPC", "Neurons", "Astrocytes",
	"Oligodendrocytes", "Microglia", "Endothelial"))

## #split by age cat and region					
tIndexes <- splitit(group)

tstatList <- lapply(tIndexes, function(i) {
	x <- rep(0, ncol(yExprs_Merge))
	x[i] <- 1
	return(genefilter::rowttests(yExprs_Merge, factor(x)))
})

numProbes=25
probeList <- lapply(tstatList, function(x) {
	y <- x[which(x[, "p.value"] < 1e-15), ]
	yUp <- y[order(y[, "dm"], decreasing = FALSE),] # signs are swapped
	rownames(yUp)[1:numProbes]
})

## filter
trainingProbes <- unique(unlist(probeList))
trainingProbes = trainingProbes[!is.na(trainingProbes)]

mergeMarkerExprs <- yExprs_Merge[trainingProbes, ]

mergeMarkerMeanExprs <- colMeans(mergeMarkerExprs)

form <- as.formula(sprintf("y ~ %s - 1", paste(levels(group),collapse = "+")))
phenoDF <- as.data.frame(model.matrix(~group - 1))
colnames(phenoDF) <- sub("^group", "", colnames(phenoDF))

## try z-score
mergeMarkerExprsZ = scale(mergeMarkerExprs)
mergeMarkerMeanExprsZ = colMeans(mergeMarkerExprsZ)

## do calibration
coefEsts <- minfi:::validationCellType(Y = mergeMarkerExprsZ, 
	pheno = phenoDF, modelFix = form)$coefEsts


## add symbol
tabOut = as.data.frame(coefEsts)
tabOut$MarkerClass =  colnames(coefEsts)[apply(coefEsts, 1, which.max)]
tabOut$Symbol = rowData(rse_geneQuake[rownames(tabOut),])$Symbol
tabOut$GencodeID = rownames(tabOut)
tabOut = tabOut[,c(13,12,11,1:10)]

## make plot
library(lattice)
theSeq = seq(-3.5,3.5,by=0.2)
mat =  coefEsts
colnames(mat) = c("iPSC", "NPC", "Fetal:Repl", "Fetal:Quies", "Adult:OPC",
	"Adult:Neuron", "Adult:Astro", "Adult:Oligo", "Adult:Microglia", "Adult:Endothelial")
rownames(mat) = rowData(rse_geneQuake[rownames(coefEsts),])$Symbol
	my.col <- colorRampPalette(c("blue","white","red"))(length(theSeq))
	
## drop 4 classes
mainIndex = which(tabOut$MarkerClass %in% c("Endothelial",
	"Fetal_quiescent", "Fetal_replicating", "iPSC", "Neurons", "NPC"))

# tabOut_eeb = tabOut[mainIndex,]
# tabOut_eeb$MarkerClass = as.factor(tabOut_eeb$MarkerClass)
# tabOut_eeb$MarkerClass = factor(tabOut_eeb$MarkerClass,levels(tabOut_eeb$MarkerClass)[c(4,6,3,2,5,1)])
# tabOut_eeb = tabOut_eeb[order(tabOut_eeb$MarkerClass),]
# mat_eeb = as.matrix(tabOut_eeb[,4:13])
# rownames(mat_eeb) = tabOut_eeb$Symbol
# pdf("singleCellGroup_exprsMatZ_main_eeb.pdf",w=36)
# print(levelplot(mat_eeb, aspect = "fill", at = theSeq,pretty=TRUE,
	# panel = panel.levelplot.raster, col.regions = my.col,
		# scales=list(x=list(rot=90,cex=2), y=list(cex=1.4)),
		# ylab = "Cell Type", xlab = ""))
# dev.off()

# write.csv(tabOut, file="singleCell_iPSC_quake_coefEsts_calibration_Zscale_geneSymbols.csv")

##################
# libd stem ######
## libd time course

load("singleCell_iPSC_quake_coefEsts_calibration_Zscale.rda")

###########
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n146.rda"))
getRPKM = recount::getRPKM

## don't use technical replicates
bioInd = which(colData(rse_gene)$Class == "Naked genomes")
rse_gene = rse_gene[,bioInd]  # n=128
colData(rse_gene)$DAY = paste0("day_",colData(rse_gene)$DAY)
colData(rse_gene)$experiment = "lieber"
rse_geneTC = rse_gene

## clean CONDITION
rse_geneTC$COND = rse_geneTC$CONDITION
rse_geneTC$COND[grep("ACC_DORSAL",rse_geneTC$COND)] = "ACC_DORSAL"
rse_geneTC$COND[grep("PLUS",rse_geneTC$COND)] = "NEURONS+ASTROS"
rse_geneTC$COND = as.factor(rse_geneTC$COND)
rse_geneTC$COND = factor(rse_geneTC$COND,levels(rse_geneTC$COND)[c(5,1,4,6,2,3)])
rse_geneTC$DAYNUM = as.numeric(ss(rse_geneTC$DAY, "_",2))

#get RPKM
yExprsTC = log2(getRPKM(rse_geneTC,length_var="Length")+1)
yExprsTC_Z = scale(yExprsTC[rownames(coefEsts),])



### heatmap of signature ###

## make plot
library(lattice)
theSeq = seq(-3.5,3.5,by=0.2)
eMat = yExprsTC_Z
colnames(eMat) = rse_geneTC$COND
rownames(eMat) = rowData(rse_geneTC[rownames(coefEsts),])$Symbol
my.col <- colorRampPalette(c("blue","white","red"))(length(theSeq))

	########################
## add symbol to blank ones (per ensembl.org)
# ensemblID
# ENSG00000254277	AC009446.1
# ENSG00000249152	LNCPRESS2
# ENSG00000280707	AL353747.4
# ENSG00000254339	AC064802.1
# ENSG00000253507	AC104257.1
# ENSG00000280511	AL591030.1
# ENSG00000256637	LINC01965
# ENSG00000278543	AC010332.2
# ENSG00000201774	RF00019
# ENSG00000253638	AF186190.1
# ENSG00000280234	AC124303.2
syms = c("AC009446.1","LNCPRESS2","AL353747.4","AC064802.1","AC104257.1",
			"AL591030.1","LINC01965","AC010332.2","RF00019","AF186190.1","AC124303.2")

ooStem = order(rse_geneTC$COND, rse_geneTC$DAYNUM)
eMat_sub = eMat[mainIndex,ooStem]
eMap_sub = rowData(rse_geneTC[rownames(coefEsts),])
eMap_sub = eMap_sub[mainIndex,]
eMap_sub$Symbol[which(eMap_sub$Symbol=="")] = syms
rownames(eMat_sub) = eMap_sub$Symbol

# pdf("humanNeuron_exprsMatZ_onlyMain_eeb.pdf",w=36,h=10)
# print(levelplot(eMat_sub, aspect = "fill", at = theSeq,pretty=TRUE,
	# panel = panel.levelplot.raster, col.regions = my.col,
		# scales=list(x=list(rot=90,cex=1.85)),	
		# ylab = "", xlab = ""))
# dev.off()




## stem-cell
stem = eMat_sub
stem = as.data.frame(t(as.data.frame(stem)))
stem$cond = ss(rownames(stem), "\\.")
stem$cond = factor(stem$cond, levels(factor(stem$cond))[c(5,1,4,6,3,2)])
levels(stem$cond) = c("Renewal","Accel Dorsal","NPC","Rosette","Neurons","Neurons + Astros")
	
## single-cell
sing = mergeMarkerExprsZ[mainIndex,]
colnames(sing) = group
rownames(sing) = rownames(eMat_sub)
sing = as.data.frame(t(as.data.frame(sing)))
sing$type = ss(rownames(sing), "\\.")
sing$type = factor(sing$type, levels(factor(sing$type))[c(5,8,4,3,10,7,1,9,6,2)])
levels(sing$type) = c("iPSC", "NPC", "Fetal:Repl", "Fetal:Quies", "Adult:OPC",
	"Adult:Neuron", "Adult:Astro", "Adult:Oligo", "Adult:Microglia", "Adult:Endothelial")

## 131 total
pdf("markerExprsZ_boxplots.pdf", h=5,w=10)
par(mar=c(9,6,2,2),cex.axis=1.2,cex.lab=1.35,cex.main=2)
for (i in 1:131) {
palette(brewer.pal(10,"RdYlGn"))
boxplot(sing[,i] ~ as.numeric(sing$type) ,outline=FALSE, xaxt="n",
		xlim = c(0.5,17.5), ylim = c(-1.2,4), las=3,
		ylab= "Standardized Expression (Z-scale)", main=colnames(sing)[i])
	points(sing[,i] ~ jitter(as.numeric(sing$type), amount=0.1),
		pch = 21, bg=pal[as.numeric(sing$type)], cex=1.3)
	axis(1, at=1:10, labels = FALSE)
	text(1:10, par("usr")[3]-.4,
		srt=40, adj= 1, xpd = TRUE,
		labels = levels(sing$type), cex=1.25)
	## add bar at means	
	means <- tapply(sing[,i], as.numeric(sing$type), mean)
	for (m in 1:length(means)) { 
		segments(x0=m-0.5, y0=means[m], x1=m+0.5, y1=means[m], col="red4", lwd=2) 
		points(x=0.6:9.6, y=means, col="red4", pch=18, cex=1.5)
		points(x=1.4:10.4, y=means, col="red4", pch=18, cex=1.5)
	}
	
palette(brewer.pal(9,"Spectral"))		
boxplot(stem[,i] ~ as.numeric(stem$cond), outline=FALSE, yaxt="n", xaxt="n",
	at = 12:17, add=TRUE, las=3)	
	points(stem[,i] ~ jitter(as.numeric(stem$cond)+11, amount=0.1), 
		pch = 21, bg=as.numeric(stem$cond)+3, cex=1.3)
	axis(1, at=12:17, labels = FALSE)
	text(12.1:17.1, par("usr")[3]-.4,
		srt=40, adj= 1, xpd = TRUE,
		labels = levels(stem$cond), cex=1.25)
	## add bar at means	
	means <- tapply(stem[,i], as.numeric(stem$cond), mean)
	for (m in 1:length(means)) { 
		segments(x0=m+10.5, y0=means[m], x1=m+11.5, y1=means[m], col="red4", lwd=2) 
		points(x=11.6:16.6, y=means, col="red4", pch=18, cex=1.5)
		points(x=12.4:17.4, y=means, col="red4", pch=18, cex=1.5)
	}
				
abline(v=11,lty=2)

}	

dev.off()





	
	