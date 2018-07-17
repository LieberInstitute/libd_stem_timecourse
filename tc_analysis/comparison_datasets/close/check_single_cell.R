###

library(SummarizedExperiment)
library(recount)
library(jaffelab)
library(RColorBrewer)

# load single cell data
load("/dcl01/ajaffe/data/lab/singleCell/close_interneurons/rse_gene_Close_interneurons_human_n1773.Rdata")

## add pheno
pheno = read.delim("/dcl01/ajaffe/data/lab/singleCell/close_interneurons/close_SraRunTable.txt",as.is=TRUE)
pheno = pheno[match(colnames(rse_gene), pheno$Run_s),]
colData(rse_gene) = cbind(colData(rse_gene), pheno)
rse_gene$DIV = as.numeric(gsub("D","", rse_gene$days_in_culture_s))
rse_gene_close = rse_gene
rse_gene_close$experiment = "close"
rse_gene_close$Type = ifelse(rse_gene_close$source_name_s == "cultured embryonic stem cells",
	"Bulk", "SingleCell")
# rse_gene_close$SampleLabel = paste0(rse_gene_close$cre_line_s, "_", rse_gene_close$DIV)
rse_gene_close$SampleLabel = paste0(rse_gene_close$cre_line_s, ": day ",rse_gene_close$DIV)
rse_gene_close$SampleLabel = as.factor(rse_gene_close$SampleLabel)
rse_gene_close$SampleLabel = factor(rse_gene_close$SampleLabel, levels = levels(rse_gene_close$SampleLabel)[c(5:8,1:4)])
	
## libd time course
MAINDIR = "/dcl01/lieber/ajaffe/lab/libd_stem_timecourse"
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n146.rda"))
## don't use technical replicates
bioInd = which(colData(rse_gene)$Class == "Naked genomes")
rse_gene = rse_gene[,bioInd]  # n=128
colData(rse_gene) = colData(rse_gene)[,c(24,26,23,15,17,20,6)]
colData(rse_gene)$DAY = paste0("day_",colData(rse_gene)$DAY)
colData(rse_gene)$experiment = "lieber"
rse_geneTC = rse_gene

## put conditions in order
rse_geneTC$CONDITION[grep("ACC_DORS",rse_geneTC$CONDITION)] = "ACC_DORSAL"
rse_geneTC$CONDITION[grep("PLUS",rse_geneTC$CONDITION)] = "NEURONS+ASTROS"
rse_geneTC$CONDITION = as.factor(rse_geneTC$CONDITION)
rse_geneTC$CONDITION = factor(rse_geneTC$CONDITION, levels = levels(rse_geneTC$CONDITION)[c(5,1,4,6,2,3)])

## rpkm function
getRPKM = recount::getRPKM

##################
## separate PCA ##
##################

###################
## based on LIBD
geneRpkmTC = getRPKM(rse_geneTC, "Length")
yRpkmTC = log2(geneRpkmTC+1)
pcaTC = prcomp(t(yRpkmTC))
pcaTC_Vars = getPcaVars(pcaTC)

# project
geneRpkmClose = getRPKM(rse_gene_close,length_var="Length")
yExprsClose_Scaled = scale(t(log2(geneRpkmClose+1)), pcaTC$center, pcaTC$scale) 
genePCs_Close_projected = yExprsClose_Scaled %*% pcaTC$rotation 

## panel A
pdf("close_projection_PCsbasedOnLIBD_A.pdf", h=6,w=7)
par(mar=c(5,6,4,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(9,"Spectral"))
plot(pcaTC$x, bg = as.numeric(factor(rse_geneTC$CONDITION))+3,pch=21,cex=2,
	xlab=paste0("PC1: ", pcaTC_Vars[1], "% Var Expl"),
    ylab=paste0("PC2: ", pcaTC_Vars[2], "% Var Expl"))
legend("bottom", paste0(levels(factor(rse_geneTC$CONDITION))), col=4:9,
       pch = 16, cex=.75)

palette(brewer.pal(4,"Oranges"))
points(genePCs_Close_projected[,1:2], pch = 22, cex=2,
	bg=as.numeric(factor(rse_gene_close$DIV)))
legend("bottomright", paste0("Day ",levels(factor(rse_gene_close$DIV))), col=1:4,
   pch = 15, cex=.85)
dev.off()



  
## panels B and C
pdf("close_projection_PCsbasedOnLIBD_B_C.pdf", h=6,w=7)
par(mar=c(8,6,4,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(4,"Oranges"))
boxplot(genePCs_Close_projected[,1] ~ rse_gene_close$DIV, outline=FALSE, ylim=c(-25,100))
points(genePCs_Close_projected[,1] ~ jitter(as.numeric(factor(rse_gene_close$DIV))),
        pch=21, 
        bg=as.numeric(factor(rse_gene_close$DIV)), cex=2)

boxplot(genePCs_Close_projected[,1] ~ rse_gene_close$SampleLabel, las=3, outline=FALSE, ylim=c(-25,100), xaxt="n")
points(genePCs_Close_projected[,1] ~ jitter(as.numeric(rse_gene_close$SampleLabel)),
        pch=21, 
        bg=as.numeric(factor(rse_gene_close$DIV)), cex=2)
axis(1, at=1:8, labels = FALSE)
text(1.3:8.3, par("usr")[3]-.07,
     srt = 60, adj= 1.2, xpd = TRUE,
     labels = levels(rse_gene_close$SampleLabel), cex=1.25)
	 
boxplot(genePCs_Close_projected[,2] ~ rse_gene_close$DIV, outline=FALSE, ylim=c(-32,30))
points(genePCs_Close_projected[,2] ~ jitter(as.numeric(factor(rse_gene_close$DIV))),
        pch=21, 
        bg=as.numeric(factor(rse_gene_close$DIV)), cex=2)
		
boxplot(genePCs_Close_projected[,2] ~ rse_gene_close$SampleLabel, las=3, outline=FALSE, ylim=c(-32,30), xaxt="n")
points(genePCs_Close_projected[,2] ~ jitter(as.numeric(rse_gene_close$SampleLabel)),
        pch=21, 
        bg=as.numeric(factor(rse_gene_close$DIV)), cex=2)
axis(1, at=1:8, labels = FALSE)
text(1.3:8.3, par("usr")[3]-.07,
     srt = 60, adj= 1.2, xpd = TRUE,
     labels = levels(rse_gene_close$SampleLabel), cex=1.25)
	  
dev.off()


#########################
# based on single cell ##
#########################

yRpkmClose= log2(geneRpkmClose+1)
pcaClose = prcomp(t(yRpkmClose))
pcaClose_Vars = getPcaVars(pcaClose)

# project
yExprsTC_Scaled = scale(t(yRpkmTC), pcaClose$center, pcaClose$scale) 
genePCs_TC_projected = yExprsTC_Scaled %*% pcaClose$rotation 

## panel A
pdf("close_projection_PCsbasedOnSingleCell_A.pdf", h=6,w=7)
par(mar=c(5,6,4,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(4,"Oranges"))
plot(pcaClose$x, bg = ifelse(rse_gene_close$Type=="Bulk","gray",as.numeric(factor(rse_gene_close$DIV))),
	pch=as.numeric(factor(rse_gene_close$cre_line_s))+21, cex=1.8,
	xlab=paste0("PC1: ", pcaClose_Vars[1], "% Var Expl"),
    ylab=paste0("PC2: ", pcaClose_Vars[2], "% Var Expl"))
legend("bottomright", c(paste0("Day ",levels(factor(rse_gene_close$DIV))),paste0(levels(factor(rse_gene_close$cre_line_s))),"Bulk", "Single Cell"),
	col=c(1:4,4,4,"gray",4), pch = c(rep(15,4),22,23,8,8), cex=.8)
   
palette(brewer.pal(9,"Spectral"))
points(genePCs_TC_projected[,1:2], pch = 21, cex=1.8,
	bg=as.numeric(factor(rse_geneTC$CONDITION))+3)
legend("topleft", paste0(levels(factor(rse_geneTC$CONDITION))), col=4:9,
       pch = 16, cex=.7)

palette(brewer.pal(4,"Oranges"))	   
plot(pcaClose$x, bg = ifelse(rse_gene_close$Type=="Bulk",as.numeric(factor(rse_gene_close$DIV)),"gray"),
	pch=as.numeric(factor(rse_gene_close$cre_line_s))+21, cex=1.8,
	xlab=paste0("PC1: ", pcaClose_Vars[1], "% Var Expl"),
    ylab=paste0("PC2: ", pcaClose_Vars[2], "% Var Expl"))
legend("bottomright", c(paste0("Day ",levels(factor(rse_gene_close$DIV))),paste0(levels(factor(rse_gene_close$cre_line_s))),"Bulk", "Single Cell"),
	col=c(1:4,4,4,4,"gray"), pch = c(rep(15,4),22,23,8,8), cex=.8)
   
palette(brewer.pal(9,"Spectral"))
points(genePCs_TC_projected[,1:2], pch = 21, cex=1.8,
	bg=as.numeric(factor(rse_geneTC$CONDITION))+3)
legend("topleft", paste0(levels(factor(rse_geneTC$CONDITION))), col=4:9,
       pch = 16, cex=.7)
dev.off()

	
## panels B and C
pdf("close_projection_PCsbasedOnSingleCell_B_C.pdf", h=6,w=7)
par(mar=c(8,6,4,2),cex.axis=1.7,cex.lab=2,cex.main=2)
palette(brewer.pal(9,"Spectral"))
boxplot(genePCs_TC_projected[,1] ~ as.numeric(ss(rse_geneTC$DAY,"_",2)) ,outline=FALSE)
points(genePCs_TC_projected[,1] ~ jitter(as.numeric(as.factor(as.numeric(ss(rse_geneTC$DAY,"_",2))))),
        pch=21, xlab = "Timecourse DIV",
        bg=as.numeric(factor(rse_geneTC$CONDITION))+3, cex=2)
legend("bottomleft", paste0(levels(factor(rse_geneTC$CONDITION))), col=4:9,
       pch = 16, cex=.8)
	   
boxplot(genePCs_TC_projected[,2] ~ as.numeric(ss(rse_geneTC$DAY,"_",2)), outline=FALSE)
points(genePCs_TC_projected[,2] ~ jitter(as.numeric(as.factor(as.numeric(ss(rse_geneTC$DAY,"_",2))))),
        pch=21, 
        bg=as.numeric(factor(rse_geneTC$CONDITION))+3, cex=2)
		
dev.off()




pdf("close_projection_PCsbasedOnSingleCell_B_C.pdf", h=6,w=10)
par(mar=c(10,6,2,2),cex.axis=1.5,cex.lab=2,cex.main=2)
## PC1
palette(brewer.pal(9,"Spectral"))
boxplot(genePCs_TC_projected[,1] ~ as.numeric(rse_geneTC$CONDITION) ,outline=FALSE, xaxt="n",
		xlim = c(0.5,13.5), ylim = c(-100,150), las=3,
		ylab= paste0("PC1: ", pcaClose_Vars[1], "% Var Expl"))
	points(genePCs_TC_projected[,1] ~ jitter(as.numeric(rse_geneTC$CONDITION)),
		pch = 21, bg=as.numeric(rse_geneTC$CONDITION)+3, cex=2)
	axis(1, at=1:6, labels = FALSE)
	text(1.1:6.1, par("usr")[3]-12,
		srt=40, adj= 1, xpd = TRUE,
		labels = levels(rse_geneTC$CONDITION), cex=1.25)

palette(brewer.pal(4,"Oranges"))	
boxplot(pcaClose$x[,1][which(rse_gene_close$Type=="Bulk")] ~ rse_gene_close$SampleLabel[which(rse_gene_close$Type=="Bulk")], outline=FALSE, yaxt="n", xaxt="n",
	at = 8:15, add=TRUE, las=3)	
	points(pcaClose$x[,1][which(rse_gene_close$Type=="Bulk")] ~ jitter(as.numeric(factor(rse_gene_close$SampleLabel))[which(rse_gene_close$Type=="Bulk")]+7,amount=0.15), pch = 22, 
		bg=as.numeric(factor(rse_gene_close$DIV))[which(rse_gene_close$Type=="Bulk")], cex=2)
	axis(1, at=8:13, labels = FALSE)
	text(8.1:13.1, par("usr")[3]-12,
		srt=40, adj= 1, xpd = TRUE,
		labels = levels(rse_gene_close$SampleLabel[which(rse_gene_close$Type=="Bulk")])[1:6], cex=1.25)
abline(v=7,lty=2)

## PC2
palette(brewer.pal(9,"Spectral"))
boxplot(genePCs_TC_projected[,2] ~ as.numeric(rse_geneTC$CONDITION) ,outline=FALSE, xaxt="n",
		xlim = c(0.5,13.5), ylim = c(-100,150), las=3,
		ylab= paste0("PC2: ", pcaClose_Vars[2], "% Var Expl"))
	points(genePCs_TC_projected[,2] ~ jitter(as.numeric(rse_geneTC$CONDITION)),
		pch = 21, bg=as.numeric(rse_geneTC$CONDITION)+3, cex=2)
	axis(1, at=1:6, labels = FALSE)
	text(1.1:6.1, par("usr")[3]-12,
		srt=40, adj= 1, xpd = TRUE,
		labels = levels(rse_geneTC$CONDITION), cex=1.25)

palette(brewer.pal(4,"Oranges"))	
boxplot(pcaClose$x[,2][which(rse_gene_close$Type=="Bulk")] ~ rse_gene_close$SampleLabel[which(rse_gene_close$Type=="Bulk")], outline=FALSE, yaxt="n", xaxt="n",
	at = 8:15, add=TRUE, las=3)	
	points(pcaClose$x[,2][which(rse_gene_close$Type=="Bulk")] ~ jitter(as.numeric(factor(rse_gene_close$SampleLabel))[which(rse_gene_close$Type=="Bulk")]+7,amount=0.15), pch = 22, 
		bg=as.numeric(factor(rse_gene_close$DIV))[which(rse_gene_close$Type=="Bulk")], cex=2)
	axis(1, at=8:13, labels = FALSE)
	text(8.1:13.1, par("usr")[3]-12,
		srt=40, adj= 1, xpd = TRUE,
		labels = levels(rse_gene_close$SampleLabel[which(rse_gene_close$Type=="Bulk")])[1:6], cex=1.25)
abline(v=7,lty=2)


dev.off()





pdf("close_projection_PCsbasedOnSingleCell_sing_B_C.pdf", h=6,w=10)
par(mar=c(10,6,2,2),cex.axis=1.5,cex.lab=2,cex.main=2)
## PC1
palette(brewer.pal(9,"Spectral"))
boxplot(genePCs_TC_projected[,1] ~ as.numeric(rse_geneTC$CONDITION) ,outline=FALSE, xaxt="n",
		xlim = c(0.5,15.5), ylim = c(-100,150), las=3,
		ylab= paste0("PC1: ", pcaClose_Vars[1], "% Var Expl"))
	points(genePCs_TC_projected[,1] ~ jitter(as.numeric(rse_geneTC$CONDITION)),
		pch = 21, bg=as.numeric(rse_geneTC$CONDITION)+3, cex=2)
	axis(1, at=1:6, labels = FALSE)
	text(1.1:6.1, par("usr")[3]-12,
		srt=40, adj= 1, xpd = TRUE,
		labels = levels(rse_geneTC$CONDITION), cex=1.25)

palette(brewer.pal(4,"Oranges"))	
boxplot(pcaClose$x[,1][which(rse_gene_close$Type=="SingleCell")] ~ rse_gene_close$SampleLabel[which(rse_gene_close$Type=="SingleCell")], outline=FALSE, yaxt="n", xaxt="n",
	at = 8:15, add=TRUE, las=3)	
	points(pcaClose$x[,1][which(rse_gene_close$Type=="SingleCell")] ~ jitter(as.numeric(factor(rse_gene_close$SampleLabel))[which(rse_gene_close$Type=="SingleCell")]+7,amount=0.15), pch = 22, 
		bg=as.numeric(factor(rse_gene_close$DIV))[which(rse_gene_close$Type=="SingleCell")], cex=2)
	axis(1, at=8:15, labels = FALSE)
	text(8.1:15.1, par("usr")[3]-12,
		srt=40, adj= 1, xpd = TRUE,
		labels = levels(rse_gene_close$SampleLabel[which(rse_gene_close$Type=="SingleCell")]), cex=1.25)
abline(v=7,lty=2)

## PC2
palette(brewer.pal(9,"Spectral"))
boxplot(genePCs_TC_projected[,2] ~ as.numeric(rse_geneTC$CONDITION) ,outline=FALSE, xaxt="n",
		xlim = c(0.5,15.5), ylim = c(-100,150), las=3,
		ylab= paste0("PC2: ", pcaClose_Vars[2], "% Var Expl"))
	points(genePCs_TC_projected[,2] ~ jitter(as.numeric(rse_geneTC$CONDITION)),
		pch = 21, bg=as.numeric(rse_geneTC$CONDITION)+3, cex=2)
	axis(1, at=1:6, labels = FALSE)
	text(1.1:6.1, par("usr")[3]-12,
		srt=40, adj= 1, xpd = TRUE,
		labels = levels(rse_geneTC$CONDITION), cex=1.25)

palette(brewer.pal(4,"Oranges"))	
boxplot(pcaClose$x[,2][which(rse_gene_close$Type=="SingleCell")] ~ rse_gene_close$SampleLabel[which(rse_gene_close$Type=="SingleCell")], outline=FALSE, yaxt="n", xaxt="n",
	at = 8:15, add=TRUE, las=3)	
	points(pcaClose$x[,2][which(rse_gene_close$Type=="SingleCell")] ~ jitter(as.numeric(factor(rse_gene_close$SampleLabel))[which(rse_gene_close$Type=="SingleCell")]+7,amount=0.15), pch = 22, 
		bg=as.numeric(factor(rse_gene_close$DIV))[which(rse_gene_close$Type=="SingleCell")], cex=2)
	axis(1, at=8:15, labels = FALSE)
	text(8.1:15.1, par("usr")[3]-12,
		srt=40, adj= 1, xpd = TRUE,
		labels = levels(rse_gene_close$SampleLabel[which(rse_gene_close$Type=="SingleCell")]), cex=1.25)
abline(v=7,lty=2)


dev.off()







pdf("close_projection_PCsbasedOnSingleCell_combined_B_C.pdf", h=6,w=14, useDingbats=FALSE)
par(mar=c(8,6,2,2),cex.axis=1.5,cex.lab=2,cex.main=2)
## PC1
palette(brewer.pal(9,"Spectral"))
boxplot(genePCs_TC_projected[,1] ~ as.numeric(rse_geneTC$CONDITION) ,outline=FALSE, xaxt="n",
		xlim = c(0.5,22.5), ylim = c(-100,150), las=3,
		ylab= paste0("PC1: ", pcaClose_Vars[1], "% Var Expl"))
	points(genePCs_TC_projected[,1] ~ jitter(as.numeric(rse_geneTC$CONDITION)),
		pch = 21, bg=as.numeric(rse_geneTC$CONDITION)+3, cex=2)
	axis(1, at=1:6, labels = FALSE)
	text(1.1:6.1, par("usr")[3]-12,
		srt=40, adj= 1, xpd = TRUE,
		labels = levels(rse_geneTC$CONDITION), cex=1.25)

palette(brewer.pal(4,"Oranges"))	
boxplot(pcaClose$x[,1][which(rse_gene_close$Type=="SingleCell")] ~ rse_gene_close$SampleLabel[which(rse_gene_close$Type=="SingleCell")], outline=FALSE, yaxt="n", xaxt="n",
	at = 8:15, add=TRUE, las=3)	
	points(pcaClose$x[,1][which(rse_gene_close$Type=="SingleCell")] ~ jitter(as.numeric(factor(rse_gene_close$SampleLabel))[which(rse_gene_close$Type=="SingleCell")]+7,amount=0.15), pch = 22, 
		bg=as.numeric(factor(rse_gene_close$DIV))[which(rse_gene_close$Type=="SingleCell")], cex=2)
	axis(1, at=8:15, labels = FALSE)
	text(8.1:15.1, par("usr")[3]-12,
		srt=40, adj= 1, xpd = TRUE,
		labels = levels(rse_gene_close$SampleLabel[which(rse_gene_close$Type=="SingleCell")]), cex=1.25)
abline(v=7,lty=2)
abline(v=16,lty=2)

palette(brewer.pal(4,"Oranges"))	
boxplot(pcaClose$x[,1][which(rse_gene_close$Type=="Bulk")] ~ rse_gene_close$SampleLabel[which(rse_gene_close$Type=="Bulk")], outline=FALSE, yaxt="n", xaxt="n",
	at = 17:24, add=TRUE, las=3)	
	points(pcaClose$x[,1][which(rse_gene_close$Type=="Bulk")] ~ jitter(as.numeric(factor(rse_gene_close$SampleLabel))[which(rse_gene_close$Type=="Bulk")]+16,amount=0.15), pch = 22, 
		bg=as.numeric(factor(rse_gene_close$DIV))[which(rse_gene_close$Type=="Bulk")], cex=2)
	axis(1, at=17:22, labels = FALSE)
	text(17.1:22.1, par("usr")[3]-12,
		srt=40, adj= 1, xpd = TRUE,
		labels = levels(rse_gene_close$SampleLabel[which(rse_gene_close$Type=="Bulk")])[1:6], cex=1.25)


axis(1, at=c(11,20), labels=c("Single-Cell","Bulk"), line=5.7, tick=FALSE, cex.axis=1.5, font=2 )

dev.off()













######################
## merge - IGNORE
geneRpkm = cbind(getRPKM(rse_gene_close, "Length"),getRPKM(rse_geneTC, "Length"))
yExprs = log2(geneRpkm+1) # transform


## do PCA
oo = order(rowSds(yExprs), decreasing=TRUE)[1:10000]
pca = prcomp(t(yExprs[oo,]))
pcaVars = getPcaVars(pca)

# days or cell types as colors
# bulk versus single cell as shapes
type = c(rse_gene_close$Type, rep("Bulk", ncol(rse_geneTC)))
## need to add color info
plot(pca$x, pch = 20+as.numeric(factor(type)), 
	xlab=paste0("PC1: ", pcaVars[1], "% Var Expl"),
    ylab=paste0("PC2: ", pcaVars[2], "% Var Expl"))