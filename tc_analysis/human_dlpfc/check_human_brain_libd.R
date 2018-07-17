#####

library(SummarizedExperiment)
library(recount)
library(jaffelab)
library(RColorBrewer)

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

## clean CONDITION
rse_geneTC$COND = rse_geneTC$CONDITION
rse_geneTC$COND[grep("ACC_DORSAL",rse_geneTC$COND)] = "ACC_DORSAL"
rse_geneTC$COND[grep("PLUS",rse_geneTC$COND)] = "NEURONS+ASTROS"
rse_geneTC$COND = as.factor(rse_geneTC$COND)
rse_geneTC$COND = factor(rse_geneTC$COND,levels(rse_geneTC$COND)[c(5,1,4,6,2,3)])

## rpkm function
getRPKM = recount::getRPKM

## human data
load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseGene_merged_n732.rda")
rse_geneDLPFC = rse_gene[,rse_gene$Dx == "Control" & 
	rse_gene$Race %in% c("AA","CAUC")]
rse_geneDLPFC$ageGroup = cut(rse_geneDLPFC$Age, c(-1,0,1,10,20,50,100))
levels(rse_geneDLPFC$ageGroup) = c("Fetal","Infant","Child","Teens","Adult","50+")

###################
## based on stem ##
###################
geneRpkmTC = getRPKM(rse_geneTC, "Length")
yRpkmTC = log2(geneRpkmTC+1)
pcaTC = prcomp(t(yRpkmTC))
pcaTC_Vars = getPcaVars(pcaTC)

# project
geneRpkmDLPFC = getRPKM(rse_geneDLPFC,length_var="Length")
yExprsDLPFC_Scaled = scale(t(log2(geneRpkmDLPFC+1)), pcaTC$center, pcaTC$scale) 
genePCs_DLPFC_projected = yExprsDLPFC_Scaled %*% pcaTC$rotation 

## panel A
pdf("brainseqDLPFC_projection_PCsbasedOnLIBD.pdf", h=6,w=7)
par(mar=c(5,6,4,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(9,"Spectral"))
plot(pcaTC$x, bg = as.numeric(factor(rse_geneTC$COND))+3,pch=21,cex=2,
	xlab=paste0("PC1: ", pcaTC_Vars[1], "% Var Expl"),
    ylab=paste0("PC2: ", pcaTC_Vars[2], "% Var Expl"))
legend("bottom", paste0(levels(factor(rse_geneTC$COND))), col=4:9,
       pch = 16, cex=.9)
	   
palette(brewer.pal(6,"BuPu"))
points(genePCs_DLPFC_projected[,1:2], pch = 22, cex=2,
	bg=factor(rse_geneDLPFC$ageGroup))
legend("topright", paste0(levels(factor(rse_geneDLPFC$ageGroup))), col=1:6,
       pch = 15, cex=.9)
dev.off()
	   
	   
	   
#########################
# based on brain  ##
#########################

yRpkmDLPFC= log2(geneRpkmDLPFC+1)
pcaDLPFC = prcomp(t(yRpkmDLPFC))
pcaDLPFC_Vars = getPcaVars(pcaDLPFC)

# project
yExprsTC_Scaled = scale(t(yRpkmTC), pcaDLPFC$center, pcaDLPFC$scale) 
genePCs_TC_projected = yExprsTC_Scaled %*% pcaDLPFC$rotation 

## panel A
pdf("brainseqDLPFC_projection_PCsbasedOnDLPFC.pdf", h=6,w=7)
par(mar=c(5,6,4,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(6,"BuPu"))
plot(pcaDLPFC$x, bg = factor(rse_geneDLPFC$ageGroup), 
	cex=2, pch=22, 
	xlab=paste0("PC1: ", pcaDLPFC_Vars[1], "% Var Expl"),
    ylab=paste0("PC2: ", pcaDLPFC_Vars[2], "% Var Expl"))
legend("bottom", paste0(levels(factor(rse_geneDLPFC$ageGroup))), col=1:6,
       pch = 15, cex=.8)
palette(brewer.pal(9,"Spectral"))
points(genePCs_TC_projected[,1:2], pch = 21,cex=2,
	bg=as.numeric(factor(rse_geneTC$COND))+3)
legend("bottomright", paste0(levels(factor(rse_geneTC$COND))), col=4:9,
       pch = 16, cex=.8)
## looks kind of weird
dev.off()


## pc1 alone??
pdf("brainseq_boxplot_projection_PCsbasedOnDLPFC.pdf", h=6,w=10, useDingbats=FALSE)
par(mar=c(8,6,2,2),cex.axis=1.5,cex.lab=2,cex.main=2)
# boxplot(genePCs_TC_projected[,1] ~ as.numeric(ss(rse_geneTC$DAY,"_",2)))
# rse_geneTC$COND = factor(rse_geneTC$COND, 
	# levels = c("RENEW", "ACC_DORSAL(2)", "NPC", "ROSETTE", 
		# "NEURONS_ALONE", "NEURONS_PLUS_ASTROS"))
palette(brewer.pal(9,"Spectral"))
boxplot(genePCs_TC_projected[,1] ~ rse_geneTC$COND, outline=FALSE, xaxt="n",
		xlim = c(0.5,13.5), ylim = c(-110,140), las=3,
		ylab= paste0("PC1: ", pcaDLPFC_Vars[1], "% Var Expl"))
	points(genePCs_TC_projected[,1] ~ jitter(as.numeric(rse_geneTC$COND),amount=0.15), pch = 21, 
		bg=as.numeric(factor(rse_geneTC$COND))+3, cex=2)
	axis(1, at=1:6, labels = FALSE)
	text(1.1:6.1, par("usr")[3]-12,
		srt=40, adj= 1, xpd = TRUE,
		labels = levels(rse_geneTC$COND), cex=1.25)
palette(brewer.pal(6,"BuPu"))
boxplot(pcaDLPFC$x[,1] ~ rse_geneDLPFC$ageGroup, outline=FALSE, yaxt="n", xaxt="n",
	at = 8:13, add=TRUE, las=3)	
	points(pcaDLPFC$x[,1] ~ jitter(as.numeric(rse_geneDLPFC$ageGroup)+7,amount=0.15), pch = 22, 
		bg=as.numeric(factor(rse_geneDLPFC$ageGroup)), cex=2)
	axis(1, at=8:13, labels = FALSE)
	text(8.1:13.1, par("usr")[3]-12,
		srt=40, adj= 1, xpd = TRUE,
		labels = levels(rse_geneDLPFC$ageGroup), cex=1.25)
abline(v=7,lty=2)

dev.off()





