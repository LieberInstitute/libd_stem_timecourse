###
library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(genefilter)
library(RColorBrewer)

MAINDIR = "/dcl01/lieber/ajaffe/lab/libd_stem_timecourse"


########################
### full samples

pheno = pheno = read.table("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/lin/lin_SraRunTable.txt", 
					header=TRUE, sep="\t", row.names=12)

load("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/lin/preprocessed_data/rse_gene_lin_dec5_n19.Rdata")
pd = colData(rse_gene)

stopifnot(identical(rownames(pd), rownames(pheno)))
pd = cbind(pd, pheno)

pd$label = paste0(pd$sample_group,": ",pd$tissue)
pd$label = as.factor(pd$label)
pd$label = factor(pd$label,levels(pd$label)[c(2,1,4,3)])


#get RPKM
yExprsTC = log2(recount::getRPKM(rse_gene,length_var="Length")+1)





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
gIndexes_tc = splitit(pd$label)
stemPropEsts_groupMeans = sapply(gIndexes_tc, 
	function(ii) colMeans(stemPropEsts[ii,]))
stemPropEsts_groupSEs = sapply(gIndexes_tc, 
	function(ii) apply(stemPropEsts[ii,],2,function(x) sd(x)/sqrt(length(x))))

round(100*stemPropEsts_groupMeans,3)


pdf("cell_type/lineplots/lin_cellTypeDecon.pdf",h=5,w=5)
palette(pal)
par(mar=c(10,6,2,2),cex.axis=1.2,cex.lab=1.6)

plot(stemPropEsts_groupMeans[1,], type="p",xaxt="n",
	lty=1, pch=19,lwd=2, ylim = c(0,1),col=1,xlab="",
	ylab = "Class Proportion",cex=1.3)
axis(1,at=seq(along=gIndexes_tc), rep("",4), las=3,cex.axis=.9)
text(1:4, par("usr")[3]-0.06, 
     srt = 55, adj= 1, xpd = TRUE, labels = levels(pd$label), cex=.9)
segments(x0=seq(along=gIndexes_tc), x1=seq(along=gIndexes_tc), col=1,lwd=2, 
	y0=stemPropEsts_groupMeans[1,] - 2*stemPropEsts_groupSEs[1,], 
	y1=stemPropEsts_groupMeans[1,] + 2*stemPropEsts_groupSEs[1,])

for(gen in seq(1,3,by=2) ) {
	start=gen
	end=gen+1
	for(i in 1:10) {
		lines(x=start:end, y=stemPropEsts_groupMeans[i,start:end], type="b",
			lty=1, pch=19,lwd=2, col=i,cex=1.3)
		segments(x0=start:end, x1=start:end, col=i,lwd=2, 
		y0=stemPropEsts_groupMeans[i,start:end] - 2*stemPropEsts_groupSEs[i,start:end], 
		y1=stemPropEsts_groupMeans[i,start:end] + 2*stemPropEsts_groupSEs[i,start:end])
	}
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
gIndexes_tc = splitit(pd$label)
stemPropEsts_groupMeans = sapply(gIndexes_tc, 
	function(ii) colMeans(stemPropEsts[ii,]))
stemPropEsts_groupSEs = sapply(gIndexes_tc, 
	function(ii) apply(stemPropEsts[ii,],2,function(x) sd(x)/sqrt(length(x))))

round(100*stemPropEsts_groupMeans,3)


pdf("brain_stage/lineplots/lin_brainStageDecon.pdf",h=5,w=5)
palette(pal)
par(mar=c(10,6,2,2),cex.axis=1.2,cex.lab=1.6)

plot(stemPropEsts_groupMeans[1,], type="p",xaxt="n",
	lty=1, pch=19,lwd=2, ylim = c(0,1),col=1,xlab="",
	ylab = "Class Proportion",cex=1.3)
axis(1,at=seq(along=gIndexes_tc), rep("",4), las=3,cex.axis=.9)
text(1:4, par("usr")[3]-0.06, 
     srt = 55, adj= 1, xpd = TRUE, labels = levels(pd$label), cex=.9)
segments(x0=seq(along=gIndexes_tc), x1=seq(along=gIndexes_tc), col=1,lwd=2, 
	y0=stemPropEsts_groupMeans[1,] - 2*stemPropEsts_groupSEs[1,], 
	y1=stemPropEsts_groupMeans[1,] + 2*stemPropEsts_groupSEs[1,])

for(gen in seq(1,3,by=2) ) {
	start=gen
	end=gen+1
	for(i in 1:8) {
		lines(x=start:end, y=stemPropEsts_groupMeans[i,start:end], type="b",
			lty=1, pch=19,lwd=2, col=i,cex=1.3)
		segments(x0=start:end, x1=start:end, col=i,lwd=2, 
		y0=stemPropEsts_groupMeans[i,start:end] - 2*stemPropEsts_groupSEs[i,start:end], 
		y1=stemPropEsts_groupMeans[i,start:end] + 2*stemPropEsts_groupSEs[i,start:end])
	}
}

dev.off()













# ### ANOVA by day, site, cell line ###

# gIndexes_tc = splitit(pd$label)
# stemPropEsts_groupMeans = sapply(gIndexes_tc, 
	# function(ii) colMeans(stemPropEsts[ii,]))
# stemPropEsts_groupSEs = sapply(gIndexes_tc, 
	# function(ii) apply(stemPropEsts[ii,],2,function(x) sd(x)/sqrt(length(x))))

# round(100*stemPropEsts_groupMeans,3)

# r = 100*stemPropEsts_groupMeans
# pd$npc = r["NPC",match(pd$label, colnames(r))]
# pd$fetal_rep = r["Fetal_replicating",match(pd$label, colnames(r))]
# pd$neuron = r["Neurons",match(pd$label, colnames(r))]

# pd$genotype_variation = relevel(pd$genotype_variation, ref="control")
# pd_early = pd[pd$DIV == "P30-40",]
# pd_late = pd[pd$DIV == "26-29",]

# fit_full = lm(npc ~ genotype_variation, data=pd_early)
# ## Does removing variable from model change it
# # fit_site = lm(npc ~ DIV + cell_line, data=pd_early)
# # fit_day = lm(npc ~ site + cell_line, data=pd_early)
# fit_line = lm(npc ~ genotype_variation, data=pd_early)
# anova(fit_site, fit_full)
# anova(fit_day, fit_full)
# anova(fit_line, fit_full)


# fit_full = lm(fetal_rep ~ DIV + cell_line + site, data=pd_early)
# ## Does removing variable from model change it
# fit_site = lm(fetal_rep ~ DIV + cell_line, data=pd_early)
# fit_day = lm(fetal_rep ~ site + cell_line, data=pd_early)
# fit_line = lm(fetal_rep ~ DIV + site, data=pd_early)
# anova(fit_site, fit_full)
# anova(fit_day, fit_full)
# anova(fit_line, fit_full)


# fit_full = lm(neuron ~ DIV + cell_line + site, data=pd_early)
# ## Does removing variable from model change it
# fit_site = lm(neuron ~ DIV + cell_line, data=pd_early)
# fit_day = lm(neuron ~ site + cell_line, data=pd_early)
# fit_line = lm(neuron ~ DIV + site, data=pd_early)
# anova(fit_site, fit_full)
# anova(fit_day, fit_full)
# anova(fit_line, fit_full)





