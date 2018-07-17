### load packages
library(jaffelab)
library(GenomicRanges)
library(rtracklayer)
library(SummarizedExperiment)
library(RColorBrewer)
library(rafalib)
library(edgeR)
library(limma)

MAINDIR = '/dcl01/lieber/ajaffe/lab/libd_stem_timecourse'
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n146.rda"))

geneCounts = assays(rse_gene)$counts
geneMap = as.data.frame(rowRanges(rse_gene))

## recount getRPKM version
getRPKM <- function(rse, length_var = 'bp_length', mapped_var = NULL) {
    mapped <- if(!is.null(mapped_var)) colData(rse)[, mapped_var] else colSums(assays(rse)$counts)      
    bg <- matrix(mapped, ncol = ncol(rse), nrow = nrow(rse), byrow = TRUE)   
    len <- if(!is.null(length_var)) rowData(rse)[, length_var] else width(rowRanges(rse))   
    wid <- matrix(len, nrow = nrow(rse), ncol = ncol(rse), byrow = FALSE)
    assays(rse)$counts / (wid/1000) / (bg/1e6)
}
geneRpkm = getRPKM(rse_gene, length_var="Length")

pd = colData(rse_gene)
pd$REP = as.numeric(pd$BIO_REP)
pd$REP[is.na(pd$REP)] = 1
pd$DAY = ordered(pd$DAY, levels=sort(as.numeric(levels(as.factor(pd$DAY)))) )
pd$DIV = as.numeric(as.character(pd$DAY))
pd$EXPERIMENT[pd$EXPERIMENT=="TIMECOURSE1"] = "TIMECOURSE"
pd$LINE = gsub('-','_',pd$LINE)

## don't use technical replicates
bioInd = which(pd$Class == "Naked genomes")
pd = pd[bioInd,]
geneRpkm = geneRpkm[,bioInd]
geneCounts = geneCounts[,bioInd]   # n=128

## don't use RENEW controls or NEURONS_ALONE in analyses
dropInd = which(pd$CONDITION %in% c("RENEW","NEURONS_ALONE"))
pd = pd[-dropInd,]
geneRpkm = geneRpkm[,-dropInd]
geneCounts = geneCounts[,-dropInd]   # n=106


pd$COND = pd$CONDITION
pd$COND[grep("NEURON",pd$COND)] = "NEURON"  ## rename NEURONS_PLUS_ASTROS
pd$COND[grep("ACC_DORSAL",pd$COND)] = "ACC_DORSAL"  ## rename ACC_DORSAL(2)
pd$COND = as.factor(pd$COND)
pd$COND = factor(pd$COND,levels(pd$COND)[c(1,3,4,2)])  ## put levels in order


###########################
#### filter expression ####
geneMap$meanExprs = rowMeans(geneRpkm)
gIndex = which(geneMap$meanExprs > 0.1)
geneRpkm = geneRpkm[gIndex,]
geneCounts = geneCounts[gIndex,]
geneMap = geneMap[gIndex,]
map = geneMap


########################################
## do voom analysis - contrasts ########

## using contrasts
mod = model.matrix(~0 + COND + LINE + totalAssignedGene, data=pd)
colnames( mod ) <- c(levels(pd$COND), colnames(mod)[-seq(length(levels(pd$COND)))])
cmtx <- makeContrasts( "NPC-ACC_DORSAL", "ROSETTE-NPC",'NEURON-ROSETTE', levels= mod)

###################
## find DE genes ##
# ## filter gene counts ##
# whichGene = which(rowMeans(geneCounts)>10)
# geneCounts = geneCounts[whichGene,]
# geneMap = geneMap[whichGene,]

dgeGene = DGEList(counts = geneCounts, genes = as.data.frame(geneMap))
dgeGene = calcNormFactors(dgeGene)
vGene = voom(dgeGene, mod, plot=FALSE)
fitGene = lmFit(vGene)
ebGene = ebayes(contrasts.fit(fitGene,cmtx))
ffGene = topTable(eBayes(contrasts.fit(fitGene,cmtx)), coef=1:3, n = nrow(dgeGene)) #by condition
ffGene$bonf.P.Val = p.adjust(ffGene$P.Value, "bonf")
ffGene = ffGene[rownames(ebGene$t),]

# voomStatsC = cbind(ffGene, ebGene$t, ebGene$p.value)
# names(voomStatsC)[19:21] = paste0("t_", names(voomStatsC)[19:21])
# names(voomStatsC)[22:24] = paste0("pvalue_", names(voomStatsC)[22:24])

# save(voomStatsC, file="rda/voom_factor_contrasts.rda")


## significance levels
pvalMat = as.matrix(ebGene$p.value)
qvalMat = pvalMat
qvalMat[,1:3] = p.adjust(pvalMat[,1:3],method="fdr") 
colnames(pvalMat) = paste0("p_",colnames(pvalMat))
colnames(qvalMat) = paste0("q_",colnames(qvalMat))

sig = cbind(pvalMat,qvalMat)



table(ffGene$"adj.P.Val" < 0.01)
# FALSE  TRUE
 # 5246 20220

colSums(sig<0.01)
# q_NPC-ACC_DORSAL   q_ROSETTE-NPC  q_NEURON-ROSETTE
# 9067              1994             12951


voomStatsC = cbind(ffGene, ebGene$t, pvalMat, qvalMat)
names(voomStatsC)[24:26] = paste0("t_", names(voomStatsC)[24:26])


save(voomStatsC, file="rda/voom_tc_genes.rda")


############################
## plot top voomStatsC

# Colors by Donor
pal = c(colorRampPalette(brewer.pal(10,"Paired")[1:2])(4), 
		colorRampPalette(brewer.pal(10,"Paired")[3:4])(3), 
		brewer.pal(10,"Paired")[9:10], brewer.pal(10,"Paired")[5:8])
pd$lineCol = pal[as.numeric(factor(pd$LINE))]


sigOrderMat = as.data.frame(apply(voomStatsC[,c(21,27:29)], 2, 
                                  function(x) order(x)[1:200]))
yExprs = as.matrix(log2(geneRpkm+1))

# #### csv of top results
# #### AD-NPC
# ffSig = voomStatsC[order(voomStatsC$"p_NPC-ACC_DORSAL"),]
# ffSig = ffSig[ffSig$"q_NPC-ACC_DORSAL" < 0.05,-c(11:12,15:19)]
# ffSig$AD_meanExprs = rowMeans(yExprs[ffSig$gencodeID,which(pd$CONDITION=="ACC_DORSAL(2)")] )
# ffSig$NPC_meanExprs = rowMeans(yExprs[ffSig$gencodeID,which(pd$CONDITION=="NPC")] )
# ffSig$ROSE_meanExprs = rowMeans(yExprs[ffSig$gencodeID,which(pd$CONDITION=="ROSETTE")] )
# ffSig$NEUR_meanExprs = rowMeans(yExprs[ffSig$gencodeID,which(pd$CONDITION=="NEURONS_PLUS_ASTROS")] )
# write.csv(ffSig, file="ordered_effects_contrasts_1_AD-NPC.csv", quote=FALSE)
# #### NPC-ROSETTE
# ffSig = voomStatsC[order(voomStatsC$"p_ROSETTE-NPC"),]
# ffSig = ffSig[ffSig$"q_ROSETTE-NPC" < 0.05,-c(11:12,15:19)]
# ffSig$AD_meanExprs = rowMeans(yExprs[ffSig$gencodeID,which(pd$CONDITION=="ACC_DORSAL(2)")] )
# ffSig$NPC_meanExprs = rowMeans(yExprs[ffSig$gencodeID,which(pd$CONDITION=="NPC")] )
# ffSig$ROSE_meanExprs = rowMeans(yExprs[ffSig$gencodeID,which(pd$CONDITION=="ROSETTE")] )
# ffSig$NEUR_meanExprs = rowMeans(yExprs[ffSig$gencodeID,which(pd$CONDITION=="NEURONS_PLUS_ASTROS")] )
# write.csv(ffSig, file="ordered_effects_contrasts_2_NPC-ROSETTE.csv", quote=FALSE)
# #### ROSETTE-NEURON
# ffSig = voomStatsC[order(voomStatsC$"p_NEURON-ROSETTE"),]
# ffSig = ffSig[ffSig$"q_NEURON-ROSETTE" < 0.05,-c(11:12,15:19)]
# ffSig$AD_meanExprs = rowMeans(yExprs[ffSig$gencodeID,which(pd$CONDITION=="ACC_DORSAL(2)")] )
# ffSig$NPC_meanExprs = rowMeans(yExprs[ffSig$gencodeID,which(pd$CONDITION=="NPC")] )
# ffSig$ROSE_meanExprs = rowMeans(yExprs[ffSig$gencodeID,which(pd$CONDITION=="ROSETTE")] )
# ffSig$NEUR_meanExprs = rowMeans(yExprs[ffSig$gencodeID,which(pd$CONDITION=="NEURONS_PLUS_ASTROS")] )
# write.csv(ffSig, file="ordered_effects_contrasts_3_ROSETTE-NEURON.csv", quote=FALSE)
# # #### OVERALL
# # ffSig = voomStatsC[order(voomStatsC$"P.Value"),]
# # ffSig = ffSig[ffSig$"adj.P.Val" < 0.0000005,-c(11:12,15:19)]
# # ffSig$AD_meanExprs = rowMeans(yExprs[ffSig$gencodeID,which(pd$CONDITION=="ACC_DORSAL(2)")] )
# # ffSig$NPC_meanExprs = rowMeans(yExprs[ffSig$gencodeID,which(pd$CONDITION=="NPC")] )
# # ffSig$ROSE_meanExprs = rowMeans(yExprs[ffSig$gencodeID,which(pd$CONDITION=="ROSETTE")] )
# # ffSig$NEUR_meanExprs = rowMeans(yExprs[ffSig$gencodeID,which(pd$CONDITION=="NEURONS_PLUS_ASTROS")] )
# # write.csv(ffSig, file="ordered_effects_contrasts_F_STAT.csv", quote=FALSE)



#### AD-NPC ####
ooL = sigOrderMat$"pvalue_NPC-ACC_DORSAL"
pdf("ordered_effects_contrasts_1_AD-NPC.pdf",h=6,w=8)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
for(i in ooL) {
 plot(as.numeric(pd$DAY),yExprs[i,] , xaxt="n",
	   pch = 21, bg=pd$lineCol,
       cex=2,xlab="Day",
       ylab="log2(Exprs + 1)",
       main = paste0(geneMap$Symbol[i], "\n", geneMap$gencodeID[i]) )
 axis(1, at=1:1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
 abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
 legend("top", paste0("p=",signif(voomStatsC$"pvalue_NPC-ACC_DORSAL"[i],3)), bg="white")
 if (i==ooL[1]) { legend("topright", levels(factor(pd$LINE)),
	pch = 15, col = pal, cex=1) }
}
dev.off()

#### NPC-ROSETTE ####
ooL = sigOrderMat$"pvalue_ROSETTE-NPC"
pdf("ordered_effects_contrasts_2_NPC-ROSETTE.pdf",h=6,w=8)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
for(i in ooL) {
 plot(as.numeric(pd$DAY),yExprs[i,] , xaxt="n",
	   pch = 21, bg=pd$lineCol,
       cex=2,xlab="Day",
       ylab="log2(Exprs + 1)",
       main = paste0(geneMap$Symbol[i], "\n", geneMap$gencodeID[i]) )
 axis(1, at=1:1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
 abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
 legend("top", paste0("p=",signif(voomStatsC$"pvalue_ROSETTE-NPC"[i],3)), bg="white")
 if (i==ooL[1]) { legend("bottomleft", levels(factor(pd$LINE)),
	pch = 15, col = pal, cex=1) }
}
dev.off()

#### ROSETTE-NEURON ####
ooL = sigOrderMat$"pvalue_NEURON-ROSETTE"
pdf("ordered_effects_contrasts_3_ROSETTE-NEURON.pdf",h=6,w=8)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
for(i in ooL) {
 plot(as.numeric(pd$DAY),yExprs[i,] , xaxt="n",
	   pch = 21, bg=pd$lineCol,
       cex=2,xlab="Day",
       ylab="log2(Exprs + 1)",
       main = paste0(geneMap$Symbol[i], "\n", geneMap$gencodeID[i]) )
 axis(1, at=1:1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
 abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
 legend("top", paste0("p=",signif(voomStatsC$"pvalue_NEURON-ROSETTE"[i],3)), bg="white")
 if (i==ooL[1]) { legend("bottomleft", levels(factor(pd$LINE)),
	pch = 15, col = pal, cex=1) }
}
dev.off()



#### Non-neighboring contrasts ####
voom_t_diff = voomStatsC
voom_t_diff$t_ROSETTE_AD = (voom_t_diff$"t_ROSETTE-NPC" + voom_t_diff$"t_NPC-ACC_DORSAL" )
voom_t_diff$t_NEURON_NPC = (voom_t_diff$"t_NEURON-ROSETTE" + voom_t_diff$"t_ROSETTE-NPC")
voom_t_diff$t_NEURON_AD = (voom_t_diff$"t_NEURON-ROSETTE" + voom_t_diff$"t_ROSETTE-NPC" + voom_t_diff$"t_NPC-ACC_DORSAL")
sigOrderMat2 = as.data.frame(apply(voom_t_diff[,c(15,25:27)], 2, 
                                  function(x) order(x, decreasing=TRUE)[1:200]))
#### NPC-NEURON ####								  
ooL = sigOrderMat2$"t_NEURON_NPC"
pdf("ordered_effects_contrasts_NPC-NEURON.pdf",h=6,w=8)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
for(i in ooL) {
 plot(as.numeric(pd$DAY),yExprs[i,] , xaxt="n",
	   pch = 21, bg=pd$lineCol,
       cex=2,xlab="Day",
       ylab="log2(Exprs + 1)",
       main = paste0(geneMap$Symbol[i], "\n", geneMap$gencodeID[i]) )
 axis(1, at=1:1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
 abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
 if (i==ooL[1]) { legend("topleft", levels(factor(pd$LINE)),
	pch = 15, col = pal, cex=1) }
}
dev.off()

#### AD-NEURON ####
ooL = sigOrderMat2$"t_NEURON_AD"
pdf("ordered_effects_contrasts_AD-NEURON.pdf",h=6,w=8)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
for(i in ooL) {
 plot(as.numeric(pd$DAY),yExprs[i,] , xaxt="n",
	   pch = 21, bg=pd$lineCol,
       cex=2,xlab="Day",
       ylab="log2(Exprs + 1)",
       main = paste0(geneMap$Symbol[i], "\n", geneMap$gencodeID[i]) )
 axis(1, at=1:1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
 abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
 if (i==ooL[1]) { legend("topleft", levels(factor(pd$LINE)),
	pch = 15, col = pal, cex=1) }
}
dev.off()

#### AD-ROSETTE ####
ooL = sigOrderMat2$"t_ROSETTE_AD"
pdf("ordered_effects_contrasts_AD-ROSETTE.pdf",h=6,w=8)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
for(i in ooL) {
 plot(as.numeric(pd$DAY),yExprs[i,] , xaxt="n",
	   pch = 21, bg=pd$lineCol,
       cex=2,xlab="Day",
       ylab="log2(Exprs + 1)",
       main = paste0(geneMap$Symbol[i], "\n", geneMap$gencodeID[i]) )
 axis(1, at=1:1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
 abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
 if (i==ooL[1]) { legend("topleft", levels(factor(pd$LINE)),
	pch = 15, col = pal, cex=1) }
}
dev.off()

#### overall ####
ooL = sigOrderMat2$"F"
pdf("ordered_effects_contrasts_F_STAT.pdf",h=6,w=8)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
for(i in ooL) {
 plot(as.numeric(pd$DAY),yExprs[i,] , xaxt="n",
	   pch = 21, bg=pd$lineCol,
       cex=2,xlab="Day",
       ylab="log2(Exprs + 1)",
       main = paste0(geneMap$Symbol[i], "\n", geneMap$gencodeID[i]) )
 axis(1, at=1:1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
 abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
 legend("top", paste0("p=",signif(voom_t_diff$"pvalue_NEURON-ROSETTE"[i],3)), bg="white")
 if (i==ooL[1]) { legend("topright", levels(factor(pd$LINE)),
	pch = 15, col = pal, cex=1) }
}
dev.off()





load("rda/voom_factor.rda")
load("rda/voom_factor_contrasts.rda")

pdf("voom_vs_voom_contrasts.pdf",h=10,w=10)
par(mfrow=c(2,2),mar=c(5,6,5,2),cex.axis=1,cex.lab=1.5,cex.main=1.8)
palette(brewer.pal(4,"Dark2"))
plot(voomStatsC$"t_NPC-ACC_DORSAL", voomStats$t_NPC,
	main = "NPC vs ACC_DORSAL",
	xlab = "t-statistic contrasts",
	ylab = "t-statistic factor")
abline(0,1,col="darkgrey")

plot((voomStatsC$"t_ROSETTE-NPC" + voomStatsC$"t_NPC-ACC_DORSAL"), 
	voomStats$t_ROSETTE,
	main = "ROSETTE vs ACC_DORSAL",
	xlab = "t-statistic contrasts (ROS-NPC + NPC-AD)",
	ylab = "t-statistic factor")
abline(0,1,col="darkgrey")

plot(voomStatsC$"t_ROSETTE-NPC", (voomStats$t_ROSETTE - voomStats$t_NPC),
	main = "ROSETTE vs NPC",
	xlab = "t-statistic contrasts",
	ylab = "t-statistic factor (ROS-NPC)")
abline(0,1,col="darkgrey")

plot(voomStatsC$"t_NEURON-ROSETTE", (voomStats$t_NEURON - voomStats$t_ROSETTE),
	main = "NEURON vs ROSETTE",
	xlab = "t-statistic contrasts",
	ylab = "t-statistic factor (NEU-ROS)")
abline(0,1,col="darkgrey")

plot((voomStatsC$"t_NEURON-ROSETTE" + voomStatsC$"t_ROSETTE-NPC"), (voomStats$t_NEURON - voomStats$t_NPC),
	main = "NEURON vs NPC",
	xlab = "t-statistic contrasts (NEU-ROS + ROS-NPC)",
	ylab = "t-statistic factor (NEU-NPC)")
abline(0,1,col="darkgrey")

plot((voomStatsC$"t_NEURON-ROSETTE" + voomStatsC$"t_ROSETTE-NPC" + voomStatsC$"t_NPC-ACC_DORSAL"), 
	voomStats$t_NEURON,
	main = "NEURON vs ACC_DORSAL",
	xlab = "t-statistic contrasts (NEU-ROS + ROS-NPC + NPC-AD)",
	ylab = "t-statistic factor")
abline(0,1,col="darkgrey")

dev.off()






############################## 
## run enrichment analysis ###
##############################

## significance levels
pvalMat = as.matrix(ebGene$p.value)
qvalMat = pvalMat
qvalMat[,1:3] = p.adjust(pvalMat[,1:3],method="fdr") 

identical(rownames(qvalMat), geneMap$gencodeID) ##check order

## separate up/down

qvalMat2 = ifelse(qvalMat < 0.05, TRUE, FALSE)
AD_NPCup = geneMap$EntrezID[qvalMat2[,1] & voomStatsC$"t_NPC-ACC_DORSAL">0]
AD_NPCdown = geneMap$EntrezID[qvalMat2[,1] & voomStatsC$"t_NPC-ACC_DORSAL"<0]

NPC_ROSEup = geneMap$EntrezID[qvalMat2[,2] & voomStatsC$"t_ROSETTE-NPC">0]
NPC_ROSEdown = geneMap$EntrezID[qvalMat2[,2] & voomStatsC$"t_ROSETTE-NPC"<0]

ROSE_NEUup = geneMap$EntrezID[qvalMat2[,3] & voomStatsC$"t_NEURON-ROSETTE">0]
ROSE_NEUdown = geneMap$EntrezID[qvalMat2[,3] & voomStatsC$"t_NEURON-ROSETTE"<0]


# split genes into condition lists
moduleGeneList = list(AD_NPCup,AD_NPCdown,NPC_ROSEup,NPC_ROSEdown, ROSE_NEUup,ROSE_NEUdown)
moduleGeneList = lapply(moduleGeneList, function(x) x[!is.na(x)])
names(moduleGeneList) = as.vector(outer(c("Up","Down"), c("AD_NPC","NPC_ROSE","ROSE_NEU"), paste, sep="_"))

## set universe of expressed genes
geneUniverse = as.character(geneMap$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

############################## 
## run enrichment analysis ###
##############################

## and adjusted
goBP_Adj <- compareCluster(moduleGeneList, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Hs.eg.db,
                ont = "BP", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = .1,
				readable= TRUE)
goMF_Adj <- compareCluster(moduleGeneList, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Hs.eg.db,
                ont = "MF", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = .1,
				readable= TRUE)
goCC_Adj <- compareCluster(moduleGeneList, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Hs.eg.db,
                ont = "CC", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = .1,
				readable= TRUE)
kegg_Adj <- compareCluster(moduleGeneList, fun = "enrichKEGG",
                universe = geneUniverse,  pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = .1)
save(goBP_Adj,goMF_Adj,goCC_Adj,kegg_Adj, file="cond_contrasts_enrich_updown.rda")

pdf("cond_contrasts_enrich_updown.pdf",h=6,w=18)
dotplot(goBP_Adj, includeAll="TRUE")
dotplot(goMF_Adj, includeAll="TRUE")
dotplot(goCC_Adj, includeAll="TRUE")
dotplot(kegg_Adj, includeAll="TRUE")
dev.off()

