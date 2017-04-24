### load packages
library(jaffelab)
library(GenomicRanges)
library(SummarizedExperiment)
library(RColorBrewer)
library(rafalib)

MAINDIR = '/dcl01/lieber/ajaffe/lab/libd_stem_timecourse'
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n146.rda"))

geneRpkm = getRPKM(rse_gene)
geneMap = rowData(rse_gene)

pd = colData(rse_gene)
pd$REP = as.numeric(pd$BIO_REP)
pd$REP[is.na(pd$REP)] = 1
pd$DAY = ordered(pd$DAY, levels=sort(as.numeric(levels(as.factor(pd$DAY)))) )
pd$DIV = as.numeric(as.character(pd$DAY))
pd$EXPERIMENT[pd$EXPERIMENT=="TIMECOURSE1"] = "TIMECOURSE"

## don't use technical replicates or RENEW controls in analyses
bioInd = which(pd$Class == "Naked genomes" & pd$CONDITION!="RENEW")
pd = pd[bioInd,]
geneRpkm = geneRpkm[,bioInd]

pd$COND = pd$CONDITION
pd$COND[grep("NEURON",pd$COND)] = "NEURON"  ## combine on/off astros
pd$COND = as.factor(pd$COND)
pd$COND = factor(pd$COND,levels(pd$COND)[c(1,3,4,2)])  ## put levels in order (AD, NPC, ROSE, NEU, RENEW)


###########################
#### filter expression ####
geneMap$meanExprs = rowMeans(geneRpkm)
gIndex = which(geneMap$meanExprs > 0.1)
geneRpkm = geneRpkm[gIndex,]
geneMap = geneMap[gIndex,]
map = geneMap

## merge data
yExprs = as.matrix(log2(geneRpkm+1))


####################################
### linear mixed effects ###########
####################################

####################################
# run modeling
library(lme4)
library(lmerTest)
library(doMC)
registerDoMC(cores=12)

pddf = as.data.frame(pd)

## using factor
options(warn=-1)
lmer_list = foreach(i = 1:nrow(yExprs)) %dopar% {
	if(i %% 5000 == 0) cat(".")
	fit = lmer(yExprs[i,] ~ COND + (1|LINE) + totalAssignedGene, data=pddf)
	list(ANOVA = anova(fit), COEF=summary(fit)$coef)
}
names(lmer_list) = rownames(yExprs)

# # fix coef?
# for(i in  which(sapply(lmer_list,function(x) ncol(x$COEF)) != 5)) {
	# x = lmer_list[[i]]$COEF
	# x = cbind(x, lmer_list[[i+1]]$COEF[,"df"],
		# pt(x[,3], lmer_list[[i+1]]$COEF[,"df"], lower=FALSE))
	# x = x[,c(1,2,4,3,5)]
	# colnames(x) = colnames( lmer_list[[i+1]]$COEF)
	# lmer_list[[i]]$COEF = x
# }
# ## fix anova?
# for(i in  which(sapply(lmer_list,function(x) ncol(x$ANOVA)) < 6)) {
	# x = lmer_list[[i]]$ANOVA
	# x = cbind(x, 54, pf(x[,4], 4, 54, lower=FALSE))
	# x = x[,c(2,3,1,5,4,6)]
	# colnames(x) = colnames( lmer_list[[i+1]]$ANOVA)
	# lmer_list[[i]]$ANOVA = x
# }
save(lmer_list, file="rda/lmer_factor.rda")


####################################
### analysis ###########

load("rda/lmer_factor.rda")


#### anova
anovaList = lapply(lmer_list, function(x) x$ANOVA)
anovaList = mclapply(anovaList, function(x) as.matrix(x)[1,5:6],mc.cores=12)
anovaMat= do.call("rbind", anovaList)

outStats = as.data.frame(anovaMat)
colnames(outStats) = c("anovaF", "anovaPvalue")
outStats$anovaBonf = p.adjust(outStats$anovaPvalue,"bonf")


## linear
slopeMat = t(sapply(lmer_list, 
	function(x) x$COEF[2:4,1]))
colnames(slopeMat) = paste0("log2FC_",levels(pd$COND)[-1])
tMat = t(sapply(lmer_list, 
	function(x) x$COEF[2:4,4]))
colnames(tMat) = paste0("t_",levels(pd$COND)[-1])
pvalMat = t(sapply(lmer_list, 
	function(x) x$COEF[2:4,5]))
colnames(pvalMat) = paste0("pvalue_",levels(pd$COND)[-1])
qvalMat= matrix(, nrow = nrow(pvalMat), ncol = 3)
qvalMat[,1:3] = p.adjust(pvalMat[,1:3],method="fdr") 
colnames(qvalMat) = paste0("qvalue_",levels(pd$COND)[-1])

outStats = cbind(outStats, slopeMat,tMat,pvalMat,qvalMat)
#outStats$logp_Condition = -log10(outStats$pvalue_Condition)


## add annotation
out = cbind(geneMap, outStats)
save(out, file="rda/lmer_outStats.rda")




####################################
### compare with voom  #############
####################################

####################################

load("rda/lmer_outStats.rda")
load("rda/voom_factor.rda")

lmerStats = out

## lmerStats vs voomStats

pdf("voom_vs_lme.pdf",h=10,w=10)
par(mfrow=c(2,2),mar=c(5,6,5,2),cex.axis=1,cex.lab=1.5,cex.main=1.8)
palette(brewer.pal(4,"Dark2"))
plot(lmerStats$t_NPC, voomStats$t_NPC,
	main = "NPC vs ACC_DORSAL",
	xlab = "t-statistic LME",
	ylab = "t-statistic voom")
abline(0,1,col="darkgrey")
plot(lmerStats$t_ROSETTE, voomStats$t_ROSETTE,
	main = "ROSETTE vs ACC_DORSAL",
	xlab = "t-statistic LME",
	ylab = "t-statistic voom")
abline(0,1,col="darkgrey")
plot(lmerStats$t_NEURON, voomStats$t_NEURON,
	main = "NEURON vs ACC_DORSAL",
	xlab = "t-statistic LME",
	ylab = "t-statistic voom")
abline(0,1,col="darkgrey")
plot(lmerStats$anovaF, voomStats$F,
	main = "ANOVA F Stat",
	xlab = "F-statistic LME",
	ylab = "F-statistic voom")
abline(0,1,col="darkgrey")
dev.off()

