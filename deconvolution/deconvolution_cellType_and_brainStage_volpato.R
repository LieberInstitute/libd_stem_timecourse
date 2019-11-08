###
library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(genefilter)
library(RColorBrewer)

MAINDIR = "/dcl01/lieber/ajaffe/lab/libd_stem_timecourse"


########################
### read in pheno
pheno = read.csv("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/volpato/volpato_pheno_bulk.csv", 
				stringsAsFactors=FALSE)
rownames(pheno) = pheno$SAMPLE_ID


########################
### full samples
load("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/volpato/bulk/rse_gene_volpato_bulk_n57.Rdata")
pd = colData(rse_gene)
pd = cbind(pheno[pd$SAMPLE_ID,], pd)

pd$DIV = paste0("Day_",pd$day)

pd$label = paste0(pd$DIV,"_",pd$cell_line)
pd$label = as.factor(pd$label)
pd$label = factor(pd$label,levels(pd$label)[c(1,3,2,4)])

pd$label2 = paste0(pd$DIV,"_",pd$site)
pd$label2 = as.factor(pd$label2)
pd$label2 = factor(pd$label2,levels(pd$label2)[c(1,6,2,7,3,8,4,9,5,10)])

pd$label3 = paste0(pd$DIV,"_",pd$site,"_",pd$cell_line)
pd$label3 = as.factor(pd$label3)
pd$label3 = factor(pd$label3,levels(pd$label3)[c(1,11,2,12,3,13,4,14,5,15,6,16,7,17,8,18,9,19,10,20)])

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


pdf("cell_type/lineplots/volpato_cellTypeDecon.pdf",h=5,w=5)
palette(pal)
par(mar=c(10,6,2,2),cex.axis=1.2,cex.lab=1.6)
plot(stemPropEsts_groupMeans[1,], type="p",xaxt="n",
	lty=1, pch=19,lwd=2, ylim = c(0,1),col=1,xlab="",
	ylab = "Class Proportion",cex=1.3)
axis(1,at=seq(along=gIndexes_tc),levels(pd$label), las=3,cex.axis=.9)
segments(x0=seq(along=gIndexes_tc), x1=seq(along=gIndexes_tc), col=1,lwd=2, 
	y0=stemPropEsts_groupMeans[1,] - 2*stemPropEsts_groupSEs[1,], 
	y1=stemPropEsts_groupMeans[1,] + 2*stemPropEsts_groupSEs[1,])
for(i in 1:10) {
	lines(stemPropEsts_groupMeans[i,1:2], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=seq(along=gIndexes_tc[1:2]), x1=seq(along=gIndexes_tc[1:2]), col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,1:2] - 2*stemPropEsts_groupSEs[i,1:2], 
	y1=stemPropEsts_groupMeans[i,1:2] + 2*stemPropEsts_groupSEs[i,1:2])
}
for(i in 1:10) {
	lines(x=3:4, y=stemPropEsts_groupMeans[i,3:4], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=3:4, x1=3:4, col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,3:4] - 2*stemPropEsts_groupSEs[i,3:4], 
	y1=stemPropEsts_groupMeans[i,3:4] + 2*stemPropEsts_groupSEs[i,3:4])
}
dev.off()

#################################### by site
## line plot
gIndexes_tc = splitit(pd$label2)
stemPropEsts_groupMeans = sapply(gIndexes_tc, 
	function(ii) colMeans(stemPropEsts[ii,]))
stemPropEsts_groupSEs = sapply(gIndexes_tc, 
	function(ii) apply(stemPropEsts[ii,],2,function(x) sd(x)/sqrt(length(x))))

round(100*stemPropEsts_groupMeans,3)


pdf("cell_type/lineplots/volpato_cellTypeDecon_bysite.pdf",h=5,w=5)
palette(pal)
par(mar=c(10,6,2,2),cex.axis=1.2,cex.lab=1.6)
plot(stemPropEsts_groupMeans[1,], type="p",xaxt="n",
	lty=1, pch=19,lwd=2, ylim = c(0,1),col=1,xlab="",
	ylab = "Class Proportion",cex=1.3)
axis(1,at=seq(along=gIndexes_tc),levels(pd$label2), las=3,cex.axis=.9)
segments(x0=seq(along=gIndexes_tc), x1=seq(along=gIndexes_tc), col=1,lwd=2, 
	y0=stemPropEsts_groupMeans[1,] - 2*stemPropEsts_groupSEs[1,], 
	y1=stemPropEsts_groupMeans[1,] + 2*stemPropEsts_groupSEs[1,])
for(i in 1:10) {
	lines(stemPropEsts_groupMeans[i,1:2], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=seq(along=gIndexes_tc[1:2]), x1=seq(along=gIndexes_tc[1:2]), col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,1:2] - 2*stemPropEsts_groupSEs[i,1:2], 
	y1=stemPropEsts_groupMeans[i,1:2] + 2*stemPropEsts_groupSEs[i,1:2])
}
for(i in 1:10) {
	lines(x=3:4, y=stemPropEsts_groupMeans[i,3:4], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=3:4, x1=3:4, col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,3:4] - 2*stemPropEsts_groupSEs[i,3:4], 
	y1=stemPropEsts_groupMeans[i,3:4] + 2*stemPropEsts_groupSEs[i,3:4])
}
for(i in 1:10) {
	lines(x=5:6, y=stemPropEsts_groupMeans[i,5:6], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=5:6, x1=5:6, col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,5:6] - 2*stemPropEsts_groupSEs[i,5:6], 
	y1=stemPropEsts_groupMeans[i,5:6] + 2*stemPropEsts_groupSEs[i,5:6])
}
for(i in 1:10) {
	lines(x=7:8, y=stemPropEsts_groupMeans[i,7:8], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=7:8, x1=7:8, col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,7:8] - 2*stemPropEsts_groupSEs[i,7:8], 
	y1=stemPropEsts_groupMeans[i,7:8] + 2*stemPropEsts_groupSEs[i,7:8])
}
for(i in 1:10) {
	lines(x=9:10, y=stemPropEsts_groupMeans[i,9:10], type="b",
		lty=1, pch=19,lwd=2, col=i,cex=1.3)
	segments(x0=9:10, x1=9:10, col=i,lwd=2, 
	y0=stemPropEsts_groupMeans[i,9:10] - 2*stemPropEsts_groupSEs[i,9:10], 
	y1=stemPropEsts_groupMeans[i,9:10] + 2*stemPropEsts_groupSEs[i,9:10])
}
dev.off()





# ###############################################################
# ############## brain stage ####################################

# pal = brewer.pal(8,"Dark2")

# load("brain_stage/iPSC_BrainSpan_coefEsts_calibration_Zscore.rda")

# yExprsTC_Z = scale(yExprsTC[rownames(coefEsts),])

# # project
# stemPropEsts = minfi:::projectCellType(yExprsTC_Z,coefEsts)
# stemPropEstsScaled = prop.table(stemPropEsts,1)

# ## line plot
# gIndexes_tc = splitit(pd$label)
# stemPropEsts_groupMeans = sapply(gIndexes_tc, 
	# function(ii) colMeans(stemPropEsts[ii,]))
# stemPropEsts_groupSEs = sapply(gIndexes_tc, 
	# function(ii) apply(stemPropEsts[ii,],2,function(x) sd(x)/sqrt(length(x))))

# round(100*stemPropEsts_groupMeans,3)


# pdf("brain_stage/lineplots/volpato_brainStageDecon.pdf",h=5,w=5)
# palette(pal)
# par(mar=c(10,6,2,2),cex.axis=1.2,cex.lab=1.6)
# plot(stemPropEsts_groupMeans[1,], type="p",xaxt="n",
	# lty=1, pch=19,lwd=2, ylim = c(0,1),col=1,xlab="",
	# ylab = "Class Proportion",cex=1.3)
# axis(1,at=seq(along=gIndexes_tc),levels(pd$label), las=3,cex.axis=.9)
# segments(x0=seq(along=gIndexes_tc), x1=seq(along=gIndexes_tc), col=1,lwd=2, 
	# y0=stemPropEsts_groupMeans[1,] - 2*stemPropEsts_groupSEs[1,], 
	# y1=stemPropEsts_groupMeans[1,] + 2*stemPropEsts_groupSEs[1,])
# for(i in 1:8) {
	# lines(stemPropEsts_groupMeans[i,1:2], type="b",
		# lty=1, pch=19,lwd=2, col=i,cex=1.3)
	# segments(x0=seq(along=gIndexes_tc[1:2]), x1=seq(along=gIndexes_tc[1:2]), col=i,lwd=2, 
	# y0=stemPropEsts_groupMeans[i,1:2] - 2*stemPropEsts_groupSEs[i,1:2], 
	# y1=stemPropEsts_groupMeans[i,1:2] + 2*stemPropEsts_groupSEs[i,1:2])
# }
# for(i in 1:8) {
	# lines(x=3:4, y=stemPropEsts_groupMeans[i,3:4], type="b",
		# lty=1, pch=19,lwd=2, col=i,cex=1.3)
	# segments(x0=3:4, x1=3:4, col=i,lwd=2, 
	# y0=stemPropEsts_groupMeans[i,3:4] - 2*stemPropEsts_groupSEs[i,3:4], 
	# y1=stemPropEsts_groupMeans[i,3:4] + 2*stemPropEsts_groupSEs[i,3:4])
# }
# dev.off()


# #################################### by site
# ## line plot
# gIndexes_tc = splitit(pd$label2)
# stemPropEsts_groupMeans = sapply(gIndexes_tc, 
	# function(ii) colMeans(stemPropEsts[ii,]))
# stemPropEsts_groupSEs = sapply(gIndexes_tc, 
	# function(ii) apply(stemPropEsts[ii,],2,function(x) sd(x)/sqrt(length(x))))

# round(100*stemPropEsts_groupMeans,3)


# pdf("brain_stage/lineplots/volpato_brainStageDecon_bysite.pdf",h=5,w=5)
# palette(pal)
# par(mar=c(10,6,2,2),cex.axis=1.2,cex.lab=1.6)
# plot(stemPropEsts_groupMeans[1,], type="p",xaxt="n",
	# lty=1, pch=19,lwd=2, ylim = c(0,1),col=1,xlab="",
	# ylab = "Class Proportion",cex=1.3)
# axis(1,at=seq(along=gIndexes_tc),levels(pd$label2), las=3,cex.axis=.9)
# segments(x0=seq(along=gIndexes_tc), x1=seq(along=gIndexes_tc), col=1,lwd=2, 
	# y0=stemPropEsts_groupMeans[1,] - 2*stemPropEsts_groupSEs[1,], 
	# y1=stemPropEsts_groupMeans[1,] + 2*stemPropEsts_groupSEs[1,])
# for(i in 1:8) {
	# lines(stemPropEsts_groupMeans[i,1:2], type="b",
		# lty=1, pch=19,lwd=2, col=i,cex=1.3)
	# segments(x0=seq(along=gIndexes_tc[1:2]), x1=seq(along=gIndexes_tc[1:2]), col=i,lwd=2, 
	# y0=stemPropEsts_groupMeans[i,1:2] - 2*stemPropEsts_groupSEs[i,1:2], 
	# y1=stemPropEsts_groupMeans[i,1:2] + 2*stemPropEsts_groupSEs[i,1:2])
# }
# for(i in 1:8) {
	# lines(x=3:4, y=stemPropEsts_groupMeans[i,3:4], type="b",
		# lty=1, pch=19,lwd=2, col=i,cex=1.3)
	# segments(x0=3:4, x1=3:4, col=i,lwd=2, 
	# y0=stemPropEsts_groupMeans[i,3:4] - 2*stemPropEsts_groupSEs[i,3:4], 
	# y1=stemPropEsts_groupMeans[i,3:4] + 2*stemPropEsts_groupSEs[i,3:4])
# }
# for(i in 1:8) {
	# lines(x=5:6, y=stemPropEsts_groupMeans[i,5:6], type="b",
		# lty=1, pch=19,lwd=2, col=i,cex=1.3)
	# segments(x0=5:6, x1=5:6, col=i,lwd=2, 
	# y0=stemPropEsts_groupMeans[i,5:6] - 2*stemPropEsts_groupSEs[i,5:6], 
	# y1=stemPropEsts_groupMeans[i,5:6] + 2*stemPropEsts_groupSEs[i,5:6])
# }
# for(i in 1:8) {
	# lines(x=7:8, y=stemPropEsts_groupMeans[i,7:8], type="b",
		# lty=1, pch=19,lwd=2, col=i,cex=1.3)
	# segments(x0=7:8, x1=7:8, col=i,lwd=2, 
	# y0=stemPropEsts_groupMeans[i,7:8] - 2*stemPropEsts_groupSEs[i,7:8], 
	# y1=stemPropEsts_groupMeans[i,7:8] + 2*stemPropEsts_groupSEs[i,7:8])
# }
# for(i in 1:8) {
	# lines(x=9:10, y=stemPropEsts_groupMeans[i,9:10], type="b",
		# lty=1, pch=19,lwd=2, col=i,cex=1.3)
	# segments(x0=9:10, x1=9:10, col=i,lwd=2, 
	# y0=stemPropEsts_groupMeans[i,9:10] - 2*stemPropEsts_groupSEs[i,9:10], 
	# y1=stemPropEsts_groupMeans[i,9:10] + 2*stemPropEsts_groupSEs[i,9:10])
# }
# dev.off()











### ANOVA by day, site, cell line ### (for cell type)

# gIndexes_tc = splitit(pd$label3)
# stemPropEsts_groupMeans = sapply(gIndexes_tc, 
	# function(ii) colMeans(stemPropEsts[ii,]))
# stemPropEsts_groupSEs = sapply(gIndexes_tc, 
	# function(ii) apply(stemPropEsts[ii,],2,function(x) sd(x)/sqrt(length(x))))

# round(100*stemPropEsts_groupMeans,3)

# r = 100*stemPropEsts_groupMeans
# pd$ipsc = r["iPSC",match(pd$label3, colnames(r))]
# pd$npc = r["NPC",match(pd$label3, colnames(r))]
# pd$fetalrep = r["Fetal_replicating",match(pd$label3, colnames(r))]
# pd$fetalqu = r["Fetal_quiescent",match(pd$label3, colnames(r))]
# pd$opc = r["OPC",match(pd$label3, colnames(r))]
# pd$neuron = r["Neurons",match(pd$label3, colnames(r))]
# pd$astro = r["Astrocytes",match(pd$label3, colnames(r))]
# pd$oligo = r["Oligodendrocytes",match(pd$label3, colnames(r))]
# pd$glia = r["Microglia",match(pd$label3, colnames(r))]
# pd$endo = r["Endothelial",match(pd$label3, colnames(r))]



stopifnot(identical(rownames(stemPropEsts), rownames(pd)))
pd = c(pd,stemPropEsts)

celltype_vars = data.frame(matrix(, nrow=10, ncol=3))
rownames(celltype_vars) = colnames(stemPropEsts)
colnames(celltype_vars) = c("Lab","Day","Line")


co = which(colnames(pd)=="iPSC")
for (i in co:ncol(pd)) {
	fit_full = lm(pd[,i] ~ DIV + cell_line + site, data=pd)
	## Does removing variable from model change it
	fit_site = lm(pd[,i] ~ DIV + cell_line, data=pd)
	fit_day = lm(pd[,i] ~ site + cell_line, data=pd)
	fit_line = lm(pd[,i] ~ DIV + site, data=pd)
	## Save p-values of anovas
	celltype_vars[i-co+1,1] = anova(fit_site, fit_full)[2,6]
	celltype_vars[i-co+1,2] = anova(fit_day, fit_full)[2,6]
	celltype_vars[i-co+1,3] = anova(fit_line, fit_full)[2,6]
}

write.csv(celltype_vars, file="volpato_nested_anova_cell_types.csv")




# fit_full = lm(npc ~ DIV + cell_line + site, data=pd)
# ## Does removing variable from model change it
# fit_site = lm(npc ~ DIV + cell_line, data=pd)
# fit_day = lm(npc ~ site + cell_line, data=pd)
# fit_line = lm(npc ~ DIV + site, data=pd)
# anova(fit_site, fit_full)
# anova(fit_day, fit_full)
# anova(fit_line, fit_full)


# fit_full = lm(fetal_rep ~ DIV + cell_line + site, data=pd)
# ## Does removing variable from model change it
# fit_site = lm(fetal_rep ~ DIV + cell_line, data=pd)
# fit_day = lm(fetal_rep ~ site + cell_line, data=pd)
# fit_line = lm(fetal_rep ~ DIV + site, data=pd)
# anova(fit_site, fit_full)
# anova(fit_day, fit_full)
# anova(fit_line, fit_full)


# fit_full = lm(neuron ~ DIV + cell_line + site, data=pd)
# ## Does removing variable from model change it
# fit_site = lm(neuron ~ DIV + cell_line, data=pd)
# fit_day = lm(neuron ~ site + cell_line, data=pd)
# fit_line = lm(neuron ~ DIV + site, data=pd)
# anova(fit_site, fit_full)
# anova(fit_day, fit_full)
# anova(fit_line, fit_full)







### Do top PCs correlate with composition values? (cell type)

pal = c("coral","gold3","darkseagreen","deepskyblue","darkorchid")

comp = stemPropEsts

pca1 = prcomp(t(yExprsTC))
pcaVars1 = getPcaVars(pca1)

## PC 1 vs 2: Explains days / cell conditions
pdf("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/volpato/pca_log2Rpkm_PC.pdf", h=6,w=6)
par(mar=c(5,6,4,2),cex.axis=2,cex.lab=2,cex.main=2)

## Lab separates in top 2 PCs
plot(-pca1$x, pch = 21, bg=pal[as.numeric(factor(pd$site))],cex=2, main="Gene PCs",
     xlab=paste0("-PC1: ", pcaVars1[1], "% Var Expl"),
     ylab=paste0("-PC2: ", pcaVars1[2], "% Var Expl"))
legend("topleft", paste0("site ",levels(factor(pd$site))),
       pch = 15, col = pal[1:5],cex=.9)

## Plot each comp esimate for top 3 PCs
for (i in 1:ncol(comp)) {
plot(comp[,i],pca1$x[,1], pch = 21, bg=pal[as.numeric(factor(pd$site))],cex=2, main="Gene PCs",
     xlab=paste0(colnames(comp)[i]," proportion"),
     ylab=paste0("PC1: ", pcaVars1[1], "% Var Expl"))
# legend("bottomright", paste0("site ",levels(factor(pd$site))),
       # pch = 15, col = pal[1:5],cex=.9)
	 
plot(comp[,i],pca1$x[,2], pch = 21, bg=pal[as.numeric(factor(pd$site))],cex=2, main="Gene PCs",
     xlab=paste0(colnames(comp)[i]," proportion"),
     ylab=paste0("PC2: ", pcaVars1[2], "% Var Expl"))
# legend("bottomright", paste0("site ",levels(factor(pd$site))),
       # pch = 15, col = pal[1:5],cex=.9)
	   
plot(comp[,i],pca1$x[,3], pch = 21, bg=pal[as.numeric(factor(pd$site))],cex=2, main="Gene PCs",
     xlab=paste0(colnames(comp)[i]," proportion"),
     ylab=paste0("PC3: ", pcaVars1[3], "% Var Expl"))
# legend("bottomright", paste0("site ",levels(factor(pd$site))),
       # pch = 15, col = pal[1:5],cex=.9)
}

dev.off()



### Below looks the same as above
# ### Only using expressed protein_coding genes (closer to their paper)
# ## Their paper uses 13.373, this uses 14,725

# geneMap = as.data.frame(rowRanges(rse_gene))
# geneMap$meanExprs = rowMeans(recount::getRPKM(rse_gene,length_var="Length"))

# rse_gene = rse_gene[which(geneMap$gene_type == "protein_coding" & geneMap$meanExprs > 0.1),]
# geneMap = as.data.frame(rowRanges(rse_gene))
# yExprsTC = log2(recount::getRPKM(rse_gene,length_var="Length")+1)

# load("cell_type/singleCell_iPSC_quake_coefEsts_calibration_Zscale.rda")
# coefEsts = coefEsts[rownames(coefEsts) %in% rownames(geneMap) ,]	# take out filtered out genes
# yExprsTC_Z = scale(yExprsTC[rownames(coefEsts),])
# # project
# stemPropEsts = minfi:::projectCellType(yExprsTC_Z,coefEsts)


# comp = stemPropEsts

# pca1 = prcomp(t(yExprsTC))
# pcaVars1 = getPcaVars(pca1)

# ## PC 1 vs 2: Explains days / cell conditions
# pdf("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/volpato/pca_log2Rpkm_PC_14725.pdf", h=6,w=6)
# par(mar=c(5,6,4,2),cex.axis=2,cex.lab=2,cex.main=2)
# palette(brewer.pal(5,"RdYlGn"))

# ## Lab separates in top 2 PCs
# plot(pca1$x, pch = 21, bg=as.numeric(factor(pd$site)),cex=2, main="Gene PCs",
     # xlab=paste0("PC1: ", pcaVars1[1], "% Var Expl"),
     # ylab=paste0("PC2: ", pcaVars1[2], "% Var Expl"))
# legend("bottomright", paste0("site ",levels(factor(pd$site))),
       # pch = 15, col = 1:8,cex=.9)

# ## Plot each comp esimate for top 3 PCs
# for (i in 1:ncol(comp)) {
# plot(comp[,i],pca1$x[,1], pch = 21, bg=as.numeric(factor(pd$site)),cex=2, main="Gene PCs",
     # xlab=paste0(colnames(comp)[i]," proportion"),
     # ylab=paste0("PC1: ", pcaVars1[1], "% Var Expl"))
# legend("bottomright", paste0("site ",levels(factor(pd$site))),
       # pch = 15, col = 1:8,cex=.9)
	 
# plot(comp[,i],pca1$x[,2], pch = 21, bg=as.numeric(factor(pd$site)),cex=2, main="Gene PCs",
     # xlab=paste0(colnames(comp)[i]," proportion"),
     # ylab=paste0("PC2: ", pcaVars1[2], "% Var Expl"))
# legend("bottomright", paste0("site ",levels(factor(pd$site))),
       # pch = 15, col = 1:8,cex=.9)
	   
# plot(comp[,i],pca1$x[,3], pch = 21, bg=as.numeric(factor(pd$site)),cex=2, main="Gene PCs",
     # xlab=paste0(colnames(comp)[i]," proportion"),
     # ylab=paste0("PC3: ", pcaVars1[3], "% Var Expl"))
# legend("bottomright", paste0("site ",levels(factor(pd$site))),
       # pch = 15, col = 1:8,cex=.9)
# }

# dev.off()






