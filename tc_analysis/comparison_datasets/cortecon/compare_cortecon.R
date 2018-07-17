###
library(minfi)
library(jaffelab)
library(SummarizedExperiment)
library(genefilter)
library(edgeR)
library(preprocessCore)
library(WGCNA)
library(RColorBrewer)

######################
## load cortecon data ## hg38
load("/dcs01/ajaffe/LIBD_StemCell/Data/processed_rawCounts_Cortecon.rda") ## for Day, other pd
pdCort2 = pd[,1:8]
load("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/cortecon/rpkmCounts_cortecon_hg38_n24.rda")
load("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/cortecon/rawCounts_cortecon_hg38_n24.rda")
pdCort = metrics
stopifnot(identical(rownames(pdCort),pdCort2$SRR))
pdCort = cbind(pdCort,pdCort2)
geneCountsCort = geneCounts
geneRpkmCort = geneRpkm
geneMapCort = geneMap
rm(exonCounts, exonRpkm, exonMap, jCounts, jRpkm, jMap, txNumReads, txTpm, txMap, pdCort2)

#######################
## load timecourse data
MAINDIR = '/dcl01/lieber/ajaffe/lab/libd_stem_timecourse'
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n146.rda"))
geneCounts = assays(rse_gene)$counts
pd = colData(rse_gene)

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
geneMap = geneMap[geneIndex,]  ## 25466 genes
############################################################


## put Cortecon genes in same order
geneMapCort = geneMapCort[rownames(geneMap),]
geneCountsCort = geneCountsCort[rownames(geneMap),]
geneRpkmCort = geneRpkmCort[rownames(geneMap),]


##############################################
## pca w/in cortecon and then project LIBD ##
##############################################

# (don't drop RENEW in libd) # n=128
## clean CONDITION
pd$COND = pd$CONDITION
pd$COND[grep("ACC_DORSAL",pd$COND)] = "ACC_DORSAL"
pd$COND = as.factor(pd$COND)
pd$COND = factor(pd$COND,levels(pd$COND)[c(5,1,4,6,2,3)])  ## put levels in order
pd$DAY = ordered(pd$DAY, levels=sort(as.numeric(levels(as.factor(pd$DAY)))) )

pdCort$Stage = as.factor(pdCort$Stage)
pdCort$Stage = factor(pdCort$Stage,levels(pdCort$Stage)[c(4,3,1,2,5)])


pcaCort = prcomp(t(log2(geneRpkmCort+1)))
pcaCort_vars = getPcaVars(pcaCort)

geneRpkmTC_Scaled = scale(t(log2(geneRpkm+1)), pcaCort$center, pcaCort$scale) 
genePCs_TC_projected = geneRpkmTC_Scaled %*% pcaCort$rotation 

pdf("cortecon_projection_PCsbasedOnCortecon.pdf", h=6,w=7)
par(mar=c(5,6,4,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(9,"PuRd"))
plot(pcaCort$x, main="Gene PCs", cex=2, 
	pch = 22,
	bg = as.numeric(factor(pdCort$Day)),
	xlab=paste0("PC1: ", pcaCort_vars[1], "% Var Expl"),
    ylab=paste0("PC2: ", pcaCort_vars[2], "% Var Expl"))
legend("topright", paste0("Day ",levels(factor(pdCort$Day))), col=1:9,
       pch=15, cex=.9)
palette(brewer.pal(9,"Spectral"))
points(genePCs_TC_projected[,1:2], pch = 21, cex=2,
	bg = as.numeric(pd$COND)+3)
legend("bottom", paste0(levels(factor(pd$COND))), col=4:9,
       pch=16, cex=.9)

palette(brewer.pal(5,"PuRd"))	   
plot(pcaCort$x, main="Gene PCs", cex=2, 
	pch = 22,
	bg = as.numeric(factor(pdCort$Stage)),
	xlab=paste0("PC1: ", pcaCort_vars[1], "% Var Expl"),
    ylab=paste0("PC2: ", pcaCort_vars[2], "% Var Expl"))
legend("topright", paste0(levels(factor(pdCort$Stage))), col=1:9,
       pch=15, cex=.9)
palette(brewer.pal(9,"Spectral"))
points(genePCs_TC_projected[,1:2], pch = 21, cex=2,
	bg = as.numeric(pd$COND)+3)
legend("bottom", paste0(levels(factor(pd$COND))), col=4:9,
       pch=16, cex=.9)
dev.off()


##############################################
## pca w/in LIBD and then project cortecon ##

pcaTC = prcomp(t(log2(geneRpkm+1)))
pcaTC_vars = getPcaVars(pcaTC)

geneRpkmCort_Scaled = scale(t(log2(geneRpkmCort+1)), pcaTC$center, pcaTC$scale) 
genePCs_Cort_projected = geneRpkmCort_Scaled %*% pcaTC$rotation 

pdf("cortecon_projection_PCsbasedOnLIBD.pdf", h=6,w=7)
par(mar=c(5,6,4,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(9,"Spectral"))
plot(pcaTC$x, main="Gene PCs", cex=2, 
	pch = 21,
	bg = as.numeric(pd$COND)+3,
	xlab=paste0("PC1: ", pcaTC_vars[1], "% Var Expl"),
    ylab=paste0("PC2: ", pcaTC_vars[2], "% Var Expl"))
legend("bottom", paste0(levels(factor(pd$COND))), col=4:9,
       pch=16, cex=.9)
palette(brewer.pal(9,"PuRd"))
points(genePCs_Cort_projected[,1:2], pch = 22, cex=2,
	bg = as.numeric(factor(pdCort$Day)))
legend("topright", paste0("Day ",levels(factor(pdCort$Day))), col=1:9,
       pch=15, cex=.9)

palette(brewer.pal(9,"Spectral"))
plot(pcaTC$x, main="Gene PCs", cex=2, 
	pch = 21,
	bg = as.numeric(pd$COND)+3,
	xlab=paste0("PC1: ", pcaTC_vars[1], "% Var Expl"),
    ylab=paste0("PC2: ", pcaTC_vars[2], "% Var Expl"))
legend("bottom", paste0(levels(factor(pd$COND))), col=4:9,
       pch=16, cex=.9)
palette(brewer.pal(5,"PuRd"))
points(genePCs_Cort_projected[,1:2], pch = 22, cex=2,
	bg = as.numeric(factor(pdCort$Stage)))
legend("topright", paste0(levels(factor(pdCort$Stage))), col=1:9,
       pch=15, cex=.9)

dev.off()


## panels B and C
pdf("cortecon_projection_PCsbasedOnLIBD_B_C.pdf", h=6,w=10)
par(mar=c(10,6,1,2),cex.axis=1.5,cex.lab=2,cex.main=2)
## PC1
palette(brewer.pal(9,"Spectral"))
boxplot(pcaTC$x[,1] ~ as.numeric(pd$COND) ,outline=FALSE, xaxt="n",
		xlim = c(0.5,16.5), ylim = c(-100,150), las=3,
		ylab= paste0("PC1: ", pcaTC_vars[1], "% Var Expl"))
	points(pcaTC$x[,1] ~ jitter(as.numeric(pd$COND)),
		pch = 21, bg=as.numeric(pd$COND)+3, cex=2)
	axis(1, at=1:6, labels = FALSE)
	text(1.1:6.1, par("usr")[3]-12,
		srt=40, adj= 1, xpd = TRUE,
		labels = levels(pd$COND), cex=1.25)
		
palette(brewer.pal(9,"PuRd"))
boxplot(genePCs_Cort_projected[,1] ~ as.factor(pdCort$Day), outline=FALSE, yaxt="n", xaxt="n",
	at = 8:16, add=TRUE, las=3)	
	points(genePCs_Cort_projected[,1] ~ jitter(as.numeric(factor(pdCort$Day))+7,amount=0.15), pch = 22, 
		bg=as.numeric(factor(pdCort$Day)), cex=2)
	axis(1, at=8:16, labels = FALSE)
	text(8.1:16.1, par("usr")[3]-15,
		srt=40, adj= 1, xpd = TRUE,
		labels = paste0("Day ",unique(pdCort$Day)), cex=1.25)
abline(v=7,lty=2)


## PC2
palette(brewer.pal(9,"Spectral"))
boxplot(pcaTC$x[,2] ~ as.numeric(pd$COND) ,outline=FALSE, xaxt="n",
		xlim = c(0.5,16.5), ylim = c(-100,150), las=3,
		ylab= paste0("PC2: ", pcaTC_vars[2], "% Var Expl"))
	points(pcaTC$x[,2] ~ jitter(as.numeric(pd$COND)),
		pch = 21, bg=as.numeric(pd$COND)+3, cex=2)
	axis(1, at=1:6, labels = FALSE)
	text(1.1:6.1, par("usr")[3]-12,
		srt=40, adj= 1, xpd = TRUE,
		labels = levels(pd$COND), cex=1.25)
		
palette(brewer.pal(9,"PuRd"))
boxplot(genePCs_Cort_projected[,2] ~ as.factor(pdCort$Day), outline=FALSE, yaxt="n", xaxt="n",
	at = 8:16, add=TRUE, las=3)	
	points(genePCs_Cort_projected[,2] ~ jitter(as.numeric(factor(pdCort$Day))+7,amount=0.15), pch = 22, 
		bg=as.numeric(factor(pdCort$Day)), cex=2)
	axis(1, at=8:16, labels = FALSE)
	text(8.1:16.1, par("usr")[3]-15,
		srt=40, adj= 1, xpd = TRUE,
		labels = paste0("Day ",unique(pdCort$Day)), cex=1.25)
abline(v=7,lty=2)
dev.off()



### p-value for PC1: corticogenesis. PC2: Renew != NPC != neuron
pdCort$PC1 = genePCs_Cort_projected[,1]
pdCort$PC2 = genePCs_Cort_projected[,2]
## PC 1 (testing for positive association)
cor.test(pdCort$Day, pdCort$PC1, alternative="greater")
cor.test(pdCort$Day, pdCort$PC1, alternative="greater")$p.val	# 3.49e-09

## PC 2 (testing Pluri/Cortical Spec/Upper Layer Formation are different)
pdCort$numStage = as.numeric(pdCort$Stage)
fit = lm(PC2 ~ poly(numStage,2), data=pdCort)
summary(fit)													#  1.65e-08

#######################
#######  compare  #####

# a) maybe just calculating eigengenes from our WGCNA in their data for the same genes 
# and b) specific timecourse gene p-values for dynamic expression 
MAINDIR = '/dcl01/lieber/ajaffe/lab/libd_stem_timecourse'
load(file.path(MAINDIR,"tc_analysis/rda/voom_tc_genes.rda"))  # voomStatsC

###################################################
#######  voom  ###################################
pdCort$Stage[which(pdCort$Stage=="Neural Differentation")] = "Pluripotency"  ## combine first 2
pdCort$Stage = as.factor(pdCort$Stage)
pdCort$Stage = factor(pdCort$Stage,levels(pdCort$Stage)[c(4,3,1,2,5)])  ## put levels in order
levels(pdCort$Stage) <- list(Plur="Pluripotency", Cort="Cortical Specification", 
							DLform="Deep Layer Formation", ULform="Upper Layer Formation")

## using contrasts
mod = model.matrix(~0 + Stage + totalAssignedGene, data=pdCort)
colnames( mod ) <- c(levels(pdCort$Stage), colnames(mod)[-seq(length(levels(pdCort$Stage)))])
cmtx <- makeContrasts( "Cort-Plur", "DLform-Cort",'ULform-DLform', levels= mod)

###################
## find DE genes ##
dgeGene = DGEList(counts = geneCountsCort, genes = as.data.frame(geneMapCort))
dgeGene = calcNormFactors(dgeGene)
vGene = voom(dgeGene, mod, plot=FALSE)
fitGene = lmFit(vGene)
ebGene = ebayes(contrasts.fit(fitGene,cmtx))
ffGene = topTable(eBayes(contrasts.fit(fitGene,cmtx)), coef=1:3, n = nrow(dgeGene)) #by condition
ffGene$bonf.P.Val = p.adjust(ffGene$P.Value, "bonf")
ffGene = ffGene[rownames(ebGene$t),]
## significance levels
pvalMat = as.matrix(ebGene$p.value)
qvalMat = pvalMat
qvalMat[,1:3] = p.adjust(pvalMat[,1:3],method="fdr") 
colnames(pvalMat) = paste0("p_",colnames(pvalMat))
colnames(qvalMat) = paste0("q_",colnames(qvalMat))

sig = cbind(pvalMat,qvalMat)
colSums(sig<1e-5)
# q_Cort-Plur   q_DLform-Cort q_ULform-DLform
# 		906              57             75
voomStatsCort = cbind(ffGene, ebGene$t, pvalMat, qvalMat)
names(voomStatsCort)[23:25] = paste0("t_", names(voomStatsCort)[23:25])

## top results
sigOrderMat = as.data.frame(apply(voomStatsCort[,c(20,26:28)], 2, 
                                  function(x) order(x)[1:200]))
yExprs = as.matrix(log2(geneRpkmCort+1))
#### overall ####
ooL = sigOrderMat$"P.Value"
pdf("ordered_effects_contrasts_cortecon.pdf",h=6,w=8)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
for(i in ooL) {
 plot(as.numeric(as.factor(pdCort$Day)),yExprs[i,] , xaxt="n",
	   pch = 21, bg="lightblue",
       cex=2,xlab="Day",
       ylab="log2(Exprs + 1)",
       main = paste0(geneMapCort$Symbol[i], "\n", geneMapCort$gencodeID[i]) )
 axis(1, at=1:length(levels(as.factor(pdCort$Day))), labels = levels(as.factor(pdCort$Day)))
 abline(v=c(2.5,4.5,6.5), col="grey", lty=2)
 legend("top", paste0("p=",signif(voomStatsCort$"P.Value"[i],3)), bg="white")
}
dev.off()

## compare LIBD and cortecon
identical(rownames(voomStatsC),rownames(voomStatsCort))
plot(voomStatsC$"F", voomStatsCort$"F")


###################################################
#######  wgcna  ###################################

################
# plot pattern #
################

load("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/tc_analysis/wgcna/rdas/stemcell_modules_signed.rda")

geneExprsCort = log2(geneRpkmCort+1) # transform


#### by moduleEigenges function
datME = moduleEigengenes(t(geneExprsCort), netAdj$colors)$eigengenes
datME = datME[,-1] ## remove unassigned

datMEmean = list()
for (i in 1:11) {  datMEmean[[i]] = unlist(lapply(split(datME[,i], as.factor(pdCort$Day)), mean))  }
MEcort = datMEmean

pdf("cluster_rotations/wgcna_patterns_cortecon_eigengenes.pdf",h=8,w=12)
par(mfrow=c(2,3), mar=c(5,5,5,2),cex.axis=1.3,cex.lab=1.8,cex.main=2)
for (i in 1:11) {
 dat = datMEmean[[i]]
 plot(1:9, dat, cex=2.5, type="b", xaxt="n", ylim=c(-max(abs(dat))*1.5, max(abs(dat))*1.5),
		  pch=16, col="darkblue",
          ylab="Eigengene", xlab="Day", main=paste0("Cluster ",i,"  (", table(netAdj$colors)[i+1],")") )
 axis(1, at=1:length(levels(as.factor(pdCort$Day))), labels = levels(as.factor(pdCort$Day)))
 abline(v=c(2.5,4.5,6.5), col="grey", lty=2)
}
dev.off()

#### by first PC
datME = data.frame(matrix(NA, nrow=24, ncol=11))
for (k in 1:11) {  datME[,k] = pcaModules[[k]]$x[,1]  }

datMEmean = list()
for (i in 1:11) {  datMEmean[[i]] = unlist(lapply(split(datME[,i], as.factor(pdCort$Day)), mean))  }
pdf("cluster_rotations/wgcna_patterns_cortecon_1st_pc.pdf",h=8,w=12)
par(mfrow=c(2,3), mar=c(5,5,5,2),cex.axis=1.3,cex.lab=1.8,cex.main=2)
for (i in 1:11) {
 dat = datMEmean[[i]]
 plot(1:9, dat, cex=2.5, type="b", xaxt="n", ylim=c(-max(abs(dat))*1.5, max(abs(dat))*1.5),
		  pch=16, col="darkblue",
          ylab="Eigengene", xlab="Day", main=paste0("Cluster ",i,"  (", table(netAdj$colors)[i+1],")") )
 axis(1, at=1:length(levels(as.factor(pdCort$Day))), labels = levels(as.factor(pdCort$Day)))
 abline(v=c(2.5,4.5,6.5), col="grey", lty=2)
}
dev.off()



################
# Plot 50 with largest rotation for each cluster
################
moduleInds = split(1:nrow(geneRpkmCort), netAdj$colors)
names(moduleInds) = paste0("Cluster_", names(moduleInds))
moduleInds$"Cluster_0" = NULL

geneRpkmModules = lapply(moduleInds, function(x) geneRpkmCort[x,])
pcaModules = lapply(geneRpkmModules, function(x) prcomp(t(log2(x+1))) )

for (k in 1:11) {

pca_k = abs(pcaModules[[k]]$rotation[,1:2])
sigOrderMat = data.frame(PC1 = order(pca_k[,1],decreasing=TRUE)[1:200])
sigOrderMat$gene = rownames(pca_k)[sigOrderMat$PC1]
sigOrderMat$Symbol = geneMapCort[sigOrderMat$gene,"Symbol"]

maxi = min(50, table(netAdj$colors)[k+1] )

## Clusters
pdf(paste0("cluster_rotations/cluster",k,"_cortecon_topRotation.pdf"),h=6,w=8)
par(mar=c(5,6,5,2),cex.axis=1.5,cex.lab=2,cex.main=2)
for (i in 1:maxi) {
gene = sigOrderMat$gene[i]
symbol = sigOrderMat$Symbol[i]
plot(log2(geneRpkmCort[gene,]+1) ~ jitter(as.numeric(as.factor(pdCort$Day)),.5), xaxt="n",
    pch = 21, bg="lightblue",
        cex=2,xlab="Day",
        ylab="log2(Exprs+1)",
		main=paste0(symbol,"\n",gene) )
  axis(1, at=1:length(levels(as.factor(pdCort$Day))), labels = levels(as.factor(pdCort$Day)))
  abline(v=c(2.5,4.5,6.5), col="grey", lty=2)
}
dev.off()
}







#### original plots (eigengenes of LIBD timecourse clusters) #####
load("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/tc_analysis/wgcna/rdas/stemcell_modules_signed.rda")
geneExprs = log2(geneRpkm+1) # transform
mod = model.matrix(~COND + LINE + totalAssignedGene, data=pd)
geneExprsAdj = cleaningY(geneExprs, mod, P=16)


#### by moduleEigenges function
datME = moduleEigengenes(t(geneExprsAdj), netAdj$colors)$eigengenes
datME = datME[,-1] ## remove unassigned

datMEmean = list()
for (i in 1:11) {
datMEmean[[i]] = unlist(lapply(split(datME[,i], pd$DAY), mean))
}
MEtc = datMEmean

pdf("cluster_rotations/wgcna_patterns_libdtimecourse_eigengenes.pdf",h=8,w=12)
par(mfrow=c(2,3), mar=c(5,5,5,2),cex.axis=1.3,cex.lab=1.8,cex.main=2)
for (i in 1:11) {
 plot(1:9, datMEmean[[i]], cex=2.5, type="b", xaxt="n", ylim=c(-.22,.22), 
		  pch=16, col="darkblue",
          ylab="Eigengene", xlab="Day", main=paste0("Cluster ",i,"  (", table(netAdj$colors)[i+1],")") )
 axis(1, at=1:length(levels(pd$DAY)), labels = levels(pd$DAY))
 abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
}
dev.off()

#### by 1st PC
moduleInds = split(1:nrow(geneExprsAdj), netAdj$colors)
names(moduleInds) = paste0("Cluster_", names(moduleInds))
moduleInds$"Cluster_0" = NULL

exprsModules = lapply(moduleInds, function(x) geneExprsAdj[x,])
pcaModules = lapply(exprsModules, function(x) prcomp(t(x)))

pdf("cluster_rotations/wgcna_patterns_libdtimecourse_1st_pc.pdf",h=8,w=12)
par(mfrow=c(2,3), mar=c(5,5,5,2),cex.axis=1.3,cex.lab=1.8,cex.main=2)
for (i in 1:11) {
datMEmean = unlist(lapply(split(pcaModules[[i]]$x[,1], pd$DAY), mean))
 plot(1:9, datMEmean, cex=2.5, type="b", xaxt="n",  ylim=c(-max(abs(datMEmean))*1.5, max(abs(datMEmean))*1.5),
		  pch=16, col="darkblue",
          ylab="Eigengene", xlab="Day", main=paste0("Cluster ",i,"  (", table(netAdj$colors)[i+1],")") )
 axis(1, at=1:length(levels(pd$DAY)), labels = levels(pd$DAY))
 abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
}
dev.off()







### both on one plot

uniqueDays = as.factor(sort(as.numeric(unique(c(levels(as.factor(pdCort$Day)),levels(pd$DAY))))))
tcInd = which(uniqueDays %in% levels(pd$DAY))
cortInd = which(uniqueDays %in% levels(as.factor(pdCort$Day)))

pdf("wgcna_patterns_comparison.pdf",h=12,w=16)
par(mfrow=c(4,3), mar=c(4,5,4,5),cex.axis=1.2,cex.lab=1.3,cex.main=1.8)
for (i in 1:11) {
 par(mar=c(3,5,4,5), cex.axis=1.2,cex.lab=1.4, cex.main=2)
 dat = MEcort[[i]]
 plot(tcInd, MEtc[[i]], pch=16, cex=2.5, type="b", ylim=c(-.22,.22), xaxt="n", col="darkred",
          ylab="LIBD Eigengene", xlab="", main=paste0("Module ",i,"  (", table(netAdj$colors)[i+1],")") )
 axis(1, at=1:length(uniqueDays), labels = uniqueDays)
 abline(v=c(6.5,9.5,12.5), col="grey", lty=2)
 text(17,0.15, srt=-90, adj = 0, labels = "Cortecon Eigengene", xpd = TRUE, cex=1.4)
  
 par(new = T)
 plot(cortInd, dat, axes=F,xlab=NA,ylab=NA, cex=2.0, type="b", ylim=c(-max(abs(dat))*1.5, max(abs(dat))*1.5),
		  pch=15, col="darkblue")
 axis(side = 4)

 if (i==1) { legend("topleft", c("LIBD","Cortecon"), pch=c(16,15), pt.cex=c(1.6,1.4), col=c("darkred","darkblue"), cex=1.4) }

 
}
dev.off()



 