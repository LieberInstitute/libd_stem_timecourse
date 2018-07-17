
## load packages
library(jaffelab)
library(GenomicRanges)
library(rtracklayer)
library(SummarizedExperiment)
library(RColorBrewer)

MAINDIR = '/dcl01/lieber/ajaffe/lab/libd_stem_timecourse'
ERCCDIR = '/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/Ercc'

load(file.path(MAINDIR,"qc_checks/libd_stemcell_timecourse_rseERCC_n157.rda"))
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n157.rda"))

ercc = assays(rse_ercc)$tpm
gCounts = assays(rse_gene)$counts
gRpkm = getRPKM(rse_gene)
gMap = rowData(rse_gene)
pd = colData(rse_gene)

############################
## pca ########

# Colors by line
pal = c(colorRampPalette(brewer.pal(10,"Paired")[1:2])(5), 
		colorRampPalette(brewer.pal(10,"Paired")[3:4])(3), 
		brewer.pal(10,"Paired")[9:10], brewer.pal(10,"Paired")[5:8])
pd$lineCol = pal[as.numeric(factor(pd$LINE))]


pca1 = prcomp(t(log2(gRpkm+1)))
pcaVars1 = getPcaVars(pca1)

pd$dayLabel = ifelse(pd$SPECIES=="HUMAN", paste0("Day: ", pd$DAY), pd$SPECIES)
pd$dayLabel = ordered(as.factor(pd$dayLabel), levels=levels(as.factor(pd$dayLabel))[c(2,4,5,7,1,3,6,8)])
levels(pd$dayLabel)[8] = "Rat Astros"

pd$COND = as.factor(pd$CONDITION)
pd$COND = factor(pd$COND,levels(pd$COND)[c(5,1,4,6,2,3)])
levels(pd$COND) = c("Renew","Accelerated Dorsal","NPC","Rosette","Neurons","Neurons plus rat astros")

## PC 1 vs 2: Explains days / cell conditions
pdf("pca_log2Rpkm_PC1_2_day_condition.pdf", h=6,w=6)
par(mar=c(5,6,4,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(8,"Spectral"))
plot(pca1$x, pch = 21, bg=as.numeric(factor(pd$dayLabel)),cex=2, main="Gene PCs",
     xlab=paste0("PC1: ", pcaVars1[1], "% Var Expl"),
     ylab=paste0("PC2: ", pcaVars1[2], "% Var Expl"))
legend("bottom", paste0(levels(factor(pd$dayLabel))),
       pch = 15, col = 1:8,cex=.9)
plot(pca1$x, pch = 21, bg=as.numeric(factor(pd$COND)),cex=2, main="Gene PCs",
     xlab=paste0("PC1: ", pcaVars1[1], "% Var Expl"),
     ylab=paste0("PC2: ", pcaVars1[2], "% Var Expl"))
legend("bottom", paste0(levels(factor(pd$COND))),
       pch = 15, col = 1:8,cex=.9)
dev.off()

## PC 3: Explains gene assignment rate
pdf("pca_log2Rpkm_PC3_geneAssign.pdf", h=6,w=6)
par(mar=c(5,6,4,2),cex.axis=1.75,cex.lab=2,cex.main=2)
plot(pd$totalAssignedGene, pca1$x[,3],
	 pch = ifelse(pd$SPECIES=="HUMAN_RAT",23,21), 
	 bg = pd$lineCol,
	 cex=2.2, main="Gene PCs",
     ylab=paste0("PC3: ", pcaVars1[3], "% Var Expl"),
     xlab= "Gene Assignment Rate")
# legend("topleft", levels(factor(pd$LINE)),
       # pch=15, col=pal, cex=.9)
dev.off()


 
###############################################################



## Colors by control

controls = c("R16-033" , "R16-054" , "R16-073")
pd$color = ifelse(pd$RNA_NO %in% controls, "#84b74a", "#808080")  # green, gray
pd$color[pd$RNA_NO=="R16-054"] = "#E69500"  # orange
pd$color[pd$RNA_NO=="R16-073"] = "#E7298A"  # pink


###############################################################
####### compare ERCC spike-ins to expected concentrations #########
### expected ercc concentration
spikeIns = read.delim("/users/ajaffe/Lieber/Projects/RNAseq/Ribozero_Compare/ercc_actual_conc.txt",
							as.is=TRUE,row.names=2)
## match row order
spikeIns = spikeIns[match(rownames(ercc),rownames(spikeIns)),]
## mix 1
mix1conc = matrix(rep(spikeIns[,"concentration.in.Mix.1..attomoles.ul."]), 
					nc = ncol(ercc), nr = nrow(ercc), byrow=FALSE)
## ercc RMSE
SE = (ercc - 10*mix1conc)^2
pd$RSME = sqrt(colMeans(SE))
# SE = (log2(ercc+1) - log2(10*mix1conc+1))^2
# pd$RSMElog = sqrt(colMeans(SE))

### quick plot
pdf('ercc_spikein_check_mix1_n157.pdf',h=12,w=18)
mypar(4,6)
for(i in 1:ncol(ercc)) {
	plot(log2(10*spikeIns[,"concentration.in.Mix.1..attomoles.ul."]+1) ~ log2(ercc[,i]+1),
		xlab="Kallisto log2(TPM+1)", ylab="Mix 1: log2(10*Concentration+1)",
		main = colnames(ercc)[i],
		xlim = c(min(log2(ercc+1)),max(log2(ercc+1))))
	abline(0, 1, lty=2)
}
dev.off()


################################################
### subset of N=21 QC samples  (R16-033 , R16-054 , R16-073)
rnum = unique(pd$RNA_NO[pd$Class=="Internal Control"])
qcInd = which(pd$RNA_NO %in% rnum)
pdSub = pd[qcInd,]
erccSub = ercc[,qcInd]
gCountsSub = gCounts[,qcInd]

pdSub$batch = paste0(pdSub$LibraryBatch,": ",pdSub$Flowcell)
pdSub$batch = ordered(as.factor(pdSub$batch), levels=levels(as.factor(pdSub$batch))[c(2:7,1)])


### Dendrogram of ERCC spike-ins, N=21 ###
dd = dist(t(log2(erccSub+1)))
hc = hclust(dd)

pdf("dendrogram_log2tpm_ercc_n21_2.pdf",h=5,w=7)
par(mar=c(8,5,2,2))
palette(brewer.pal(8,"Dark2"))
myplclust(hc, lab.col=as.numeric(factor(pdSub$batch)), xlab="",
          lab = paste0(pdSub$RNA_NO,"_",pdSub$Library), main = "")
legend("topright", paste0(levels(factor(pdSub$batch))),
       pch=15, col=1:8, cex=.75)  	
myplclust(hc, lab.col=pdSub$color, xlab="",
          lab = paste0(pdSub$RNA_NO,"_",pdSub$Library), main = "")
legend("topright", paste0(levels(factor(pdSub$RNA_NO))),
       pch=15, col=levels(factor(pdSub$color)), cex=1)  		   
dev.off()

### Dendrogram of gene counts, N=21 ###
dd = dist(t(log2(gCountsSub+1)))
hc = hclust(dd)

pdf("dendrogram_log2_genecounts_n21.pdf",h=5,w=7)
par(mar=c(8,5,2,2))
palette(brewer.pal(8,"Dark2"))
myplclust(hc, lab.col=pdSub$color, xlab="",
          lab = paste0(pdSub$RNA_NO,"_",pdSub$Library), main = "")
legend("topright", paste0(levels(factor(pdSub$RNA_NO))),
       pch=15, col=levels(factor(pdSub$color)), cex=1)  		   
dev.off()


################################################
### 6 boxplots of mapping rate, gene assignment rate, 
### and ERCC RMSE by sequencing batch and condition, N=157

## order CONDITION categories
pd$COND = ordered(as.factor(pd$CONDITION), levels=levels(as.factor(pd$CONDITION))[c(5,1,4,6,2,3)])
lablist = c("Renew","Accelerated Dorsal","NPC","Rosette","Neurons","Neurons + Astros")

## order so colored points are plotted last/on top
pd = pd[order(pd$color),]

##### map rate #####

pdf("boxplot_mapRate_vs_batch_n157.pdf",h=4.5,w=5)
par(mar=c(6,6,2,2))
palette(brewer.pal(8,"Dark2"))
boxplot(pd$overallMapRate ~ pd$Library,
        cex.axis=1.2,cex.lab=1.5,las=3,
        ylab="Mapping Rate",outline=FALSE,
        ylim = c(0,1))
points(pd$overallMapRate ~ jitter(as.numeric(factor(pd$Library))),
       pch=as.numeric(factor(pd$SPECIES))+20, 
	   bg=pd$color,cex=1.5)
legend("bottomleft", c("HUMAN", "HUMAN+RAT"),
       pch=c(16,15), cex=.75)
dev.off()

pdf("boxplot_mapRate_vs_condition_n157.pdf",h=4.5,w=5)
par(mar=c(6,6,2,2))
palette(brewer.pal(8,"Dark2"))
boxplot(pd$overallMapRate ~ pd$COND,
        cex.axis=1.2,cex.lab=1.5,las=3,
        ylab="Mapping Rate",outline=FALSE,
        ylim = c(0,1), xaxt="n")
axis(1, at=c(1:6), labels=NA)
text(seq(1, 6, by=1), par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]), 
		labels = lablist, srt = 45, adj = 1, xpd = TRUE)
points(pd$overallMapRate ~ jitter(as.numeric(factor(pd$COND))),
       pch=as.numeric(factor(pd$SPECIES))+20, 
	   bg=pd$color, cex=1.5)
# legend("bottomleft", c("Non-control","R16-033","R16-054","R16-073"),
       # pch=15, col=levels(factor(pd$color)), cex=1)
dev.off()

##### gene assign rate #####

pdf("boxplot_geneAssignRate_vs_batch_n157.pdf",h=4.5,w=5)
par(mar=c(6,6,2,2))
palette(brewer.pal(8,"Dark2"))
boxplot(pd$totalAssignedGene ~ pd$Library,
        cex.axis=1.2,cex.lab=1.5,las=3,
        ylab="Gene Assignment Rate",outline=FALSE,
        ylim = c(0,1))
points(pd$totalAssignedGene ~ jitter(as.numeric(factor(pd$Library))),
       pch=as.numeric(factor(pd$SPECIES))+20, 
	   bg=pd$color, cex=1.5)
legend("bottomleft", c("HUMAN", "HUMAN+RAT"),
       pch=c(16,15), cex=.75)
dev.off()

pdf("boxplot_geneAssignRate_vs_condition_n157.pdf",h=4.5,w=5)
par(mar=c(6,6,2,2))
palette(brewer.pal(8,"Dark2"))
boxplot(pd$totalAssignedGene ~ pd$COND,
        cex.axis=1.2,cex.lab=1.5,las=3,
        ylab="Gene Assignment Rate",outline=FALSE,
        ylim = c(0,1), xaxt="n")
axis(1, at=c(1:6), labels=NA)
text(seq(1, 6, by=1), par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]), 
		labels = lablist, srt = 45, adj = 1, xpd = TRUE)
points(pd$totalAssignedGene ~ jitter(as.numeric(factor(pd$COND))),
       pch=as.numeric(factor(pd$SPECIES))+20, 
	   bg=pd$color, cex=1.5)
# legend("bottomleft", c("HUMAN", "HUMAN+RAT"),
       # pch=c(16,15), cex=1)
dev.off()

##### ercc RSME #####

pdf("boxplot_erccRSME_vs_batch_n157.pdf",h=4.5,w=5)
par(mar=c(6,6,2,2))
palette(brewer.pal(8,"Dark2"))
boxplot(pd$RSME ~ pd$Library,
        cex.axis=1.2,cex.lab=1.5,las=3,
        ylab="ERCC RMSE",outline=FALSE,
        ylim = c(min(pd$RSME)-500, max(pd$RSME)+700))
points(pd$RSME ~ jitter(as.numeric(factor(pd$Library))),
       pch=as.numeric(factor(pd$SPECIES))+20, 
	   bg=pd$color,cex=1.5)
dev.off()

pdf("boxplot_erccRSME_vs_condition_n157.pdf",h=4.5,w=5)
par(mar=c(6,6,2,2))
palette(brewer.pal(8,"Dark2"))
boxplot(pd$RSME ~ pd$COND,
        cex.axis=1.2,cex.lab=1.5,las=3,
        ylab="ERCC RMSE",outline=FALSE, xaxt="n",
        ylim = c(min(pd$RSME)-500, max(pd$RSME)+700))
axis(1, at=c(1:6), labels=NA)
text(seq(1, 6, by=1), par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]), 
		labels = lablist, srt = 45, adj = 1, xpd = TRUE)
points(pd$RSME ~ jitter(as.numeric(factor(pd$COND))),
       pch=as.numeric(factor(pd$SPECIES))+20, 
	   bg=pd$color,cex=1.5)
dev.off()




#############################################
### gene assignment rate without Unassigned_Multimapping

## newer metrics file with assigned and unassigned numbers
pd2 = read.csv("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/read_and_alignment_metrics_AZpilot_jan6.hg38.csv", header=TRUE, row.names=1)
## subset from 506 to 157
pd2 = pd2[pd$SAMPLE_ID,]
pd2$color = pd$color
pd2$SPECIES = pd$SPECIES

## order CONDITION categories
pd$COND = ordered(as.factor(pd$CONDITION), levels=levels(as.factor(pd$CONDITION))[c(5,1,4,6,2,3)])
lablist = c("Renew","Accelerated Dorsal","NPC","Rosette","Neurons","Neurons + Astros")
pd2$COND = pd$COND

## recalculate totalAssignedGene
pd2$assignedGene = pd2$gene_Assigned/(pd2$gene_Assigned + pd2$gene_Unassigned_Ambiguity + pd2$gene_Unassigned_NoFeatures)

## order so colored points are plotted last/on top
pd2 = pd2[order(pd2$color),]

pdf("boxplot_geneAssignRate_vs_condition_n157_without_multimapping.pdf",h=4.5,w=5)
par(mar=c(6,6,2,2))
palette(brewer.pal(8,"Dark2"))
boxplot(pd2$assignedGene ~ pd2$COND,
        cex.axis=1.2,cex.lab=1.5,las=3,
        ylab="Gene Assignment Rate",outline=FALSE,
        ylim = c(0,1), xaxt="n")
axis(1, at=c(1:6), labels=NA)
text(seq(1, 6, by=1), par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]), 
		labels = lablist, srt = 45, adj = 1, xpd = TRUE)
points(pd2$assignedGene ~ jitter(as.numeric(factor(pd2$COND))),
       pch=as.numeric(factor(pd2$SPECIES))+20, 
	   bg=pd2$color, cex=1.5)
dev.off()



# ################################################
# ### bar charts of neuron samples, different species

# library(ggplot2)
# library(reshape)
# RATDIR = "/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/ALIGN_RAT_MOUSE_n157"
# metrics = read.csv(file.path(RATDIR,"read_and_alignment_metrics_COMBINED_n157.csv"),
			# row.names=1, stringsAsFactors=FALSE)

# ## only need neuronal samples
# neuronsMap = metrics[grep("NEURONS",metrics$CONDITION),
			# c("CONDITION", "overallMapRateHUMAN","overallMapRateRAT","overallMapRateMOUSE")]	
# neuronsMito = metrics[grep("NEURONS",metrics$CONDITION),
			# c("CONDITION", "mitoRateHUMAN", "mitoRateRAT", "mitoRateMOUSE")]	
# neuronsGene = metrics[grep("NEURONS",metrics$CONDITION),
			# c("CONDITION", "totalAssignedGeneHUMAN", "totalAssignedGeneRAT", "totalAssignedGeneMOUSE")]	
# neuronsRibo = metrics[grep("NEURONS",metrics$CONDITION),
			# c("CONDITION", "rRNA_rateHUMAN", "rRNA_rateRAT", "rRNA_rateMOUSE")]				
# colnames(neuronsMap) = c("CONDITION","HUMAN\nhg38", "RAT\nrn6","MOUSE\nmm10")
# colnames(neuronsMito) = colnames(neuronsGene) = colnames(neuronsRibo) = colnames(neuronsMap)

# ##### Map Rate
# ## get means stratified by on/off astrocytes				
# means = aggregate(neuronsMap,by=list(neuronsMap$CONDITION),mean)
# means = means[,-2]
# means.long = melt(means,id.vars="Group.1")
# names(means.long) = c("Condition","genome","rate")
# ## get std devs				
# sds = aggregate(neuronsMap,by=list(neuronsMap$CONDITION),sd)
# sds = sds[,-2]
# sds.long = melt(sds,id.vars="Group.1")
# names(sds.long) = c("Condition","genome","rate")
# ## barplot by genome and condition
# pdf("neuron_samps_maprate_by_species.pdf", h=6, w=6)
# ggplot(means.long, aes(x=genome, y=rate, fill=factor(Condition)))+
  # geom_bar(stat="identity", position="dodge")+
  # geom_errorbar(aes(ymin=means.long$rate-sds.long$rate, ymax=means.long$rate+sds.long$rate),
        # width=0.1, position=position_dodge(.9))+
  # scale_fill_manual(values=c("#70b670","#3290d8"), name="Astrocytes", labels=c("Off","On"))+
  # xlab("")+
  # ylab("Alignment Rate")+
  # ylim(0,1)+
  # theme(panel.grid.major.x = element_blank(),
		# axis.text=element_text(size=12, face="bold"),
		# axis.title=element_text(size=16) )
# dev.off()


# ##### Mito rate
# ## get means stratified by on/off astrocytes				
# means = aggregate(neuronsMito,by=list(neuronsMito$CONDITION),mean)
# means = means[,-2]
# means.long = melt(means,id.vars="Group.1")
# names(means.long) = c("Condition","genome","rate")
# ## get std devs				
# sds = aggregate(neuronsMito,by=list(neuronsMito$CONDITION),sd)
# sds = sds[,-2]
# sds.long = melt(sds,id.vars="Group.1")
# names(sds.long) = c("Condition","genome","rate")
# ## barplot by genome and condition
# pdf("neuron_samps_mitorate_by_species.pdf", h=6, w=6)
# ggplot(means.long, aes(x=genome, y=rate, fill=factor(Condition)))+
  # geom_bar(stat="identity", position="dodge")+
  # geom_errorbar(aes(ymin=means.long$rate-sds.long$rate, ymax=means.long$rate+sds.long$rate),
        # width=0.1, position=position_dodge(.9))+
  # scale_fill_manual(values=c("#70b670","#3290d8"), name="Astrocytes", labels=c("Off","On"))+
  # xlab("")+
  # ylab("Mitochondrial Map Rate")+
  # ylim(0,.041)+
  # theme(panel.grid.major.x = element_blank(),
		# axis.text=element_text(size=12, face="bold"),
		# axis.title=element_text(size=16) )
# dev.off()


# ##### Gene Assignment rate
# ## get means stratified by on/off astrocytes				
# means = aggregate(neuronsGene,by=list(neuronsGene$CONDITION),mean)
# means = means[,-2]
# means.long = melt(means,id.vars="Group.1")
# names(means.long) = c("Condition","genome","rate")
# ## get std devs				
# sds = aggregate(neuronsGene,by=list(neuronsGene$CONDITION),sd)
# sds = sds[,-2]
# sds.long = melt(sds,id.vars="Group.1")
# names(sds.long) = c("Condition","genome","rate")
# ## barplot by genome and condition
# pdf("neuron_samps_geneassignrate_by_species.pdf", h=6, w=6)
# ggplot(means.long, aes(x=genome, y=rate, fill=factor(Condition)))+
  # geom_bar(stat="identity", position="dodge")+
  # # geom_errorbar(aes(ymin=means.long$rate-sds.long$rate, ymax=means.long$rate+sds.long$rate),
        # # width=0.1, position=position_dodge(.9))+
  # scale_fill_manual(values=c("#70b670","#3290d8"), name="Astrocytes", labels=c("Off","On"))+
  # xlab("")+
  # ylab("Gene Assignment Rate")+
  # ylim(0,1)+
  # theme(panel.grid.major.x = element_blank(),
		# axis.text=element_text(size=12, face="bold"),
		# axis.title=element_text(size=16) )
# dev.off()


# ##### rRNA rate
# ## get means stratified by on/off astrocytes				
# means = aggregate(neuronsRibo,by=list(neuronsRibo$CONDITION),mean)
# means = means[,-2]
# means.long = melt(means,id.vars="Group.1")
# names(means.long) = c("Condition","genome","rate")
# ## get std devs				
# sds = aggregate(neuronsRibo,by=list(neuronsRibo$CONDITION),sd)
# sds = sds[,-2]
# sds.long = melt(sds,id.vars="Group.1")
# names(sds.long) = c("Condition","genome","rate")
# ## barplot by genome and condition
# pdf("neuron_samps_rRNArate_by_species.pdf", h=6, w=6)
# ggplot(means.long, aes(x=genome, y=rate, fill=factor(Condition)))+
  # geom_bar(stat="identity", position="dodge")+
  # geom_errorbar(aes(ymin=means.long$rate-sds.long$rate, ymax=means.long$rate+sds.long$rate),
        # width=0.1, position=position_dodge(.9))+
  # scale_fill_manual(values=c("#70b670","#3290d8"), name="Astrocytes", labels=c("Off","On"))+
  # xlab("")+
  # ylab("rRNA Assignment Rate")+
  # ylim(0,.035)+
  # theme(panel.grid.major.x = element_blank(),
		# axis.text=element_text(size=12, face="bold"),
		# axis.title=element_text(size=16) )
# dev.off()




# ################################################
# ### bar charts of human map rate from rat and mouse hits

# library(ggplot2)
# library(reshape)
# RATDIR = "/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/ALIGN_RAT_MOUSE_n157/bamtofastq"
# metrics = read.csv(file.path(RATDIR,"read_and_alignment_metrics_COMBINED_n157.csv"),
			# row.names=1, stringsAsFactors=FALSE)

# ## only need neuronal samples
# neuronsMap = metrics[grep("NEURONS",metrics$CONDITION),
			# c("CONDITION","overallMapRateRAT","overallMapRateMOUSE")]	
		
# colnames(neuronsMap) = c("CONDITION","RAT\nrn6","MOUSE\nmm10")

# ##### Map Rate
# ## get means stratified by on/off astrocytes				
# means = aggregate(neuronsMap,by=list(neuronsMap$CONDITION),mean)
# means = means[,-2]
# means.long = melt(means,id.vars="Group.1")
# names(means.long) = c("Condition","genome","rate")
# ## get std devs				
# sds = aggregate(neuronsMap,by=list(neuronsMap$CONDITION),sd)
# sds = sds[,-2]
# sds.long = melt(sds,id.vars="Group.1")
# names(sds.long) = c("Condition","genome","rate")
# ## barplot by genome and condition
# pdf("neuron_samps_maprate_mouse_rat_mapped_to_human.pdf", h=6, w=6)
# ggplot(means.long, aes(x=genome, y=rate, fill=factor(Condition)))+
  # geom_bar(stat="identity", position="dodge")+
  # geom_errorbar(aes(ymin=means.long$rate-sds.long$rate, ymax=means.long$rate+sds.long$rate),
        # width=0.1, position=position_dodge(.9))+
  # scale_fill_manual(values=c("#70b670","#3290d8"), name="Astrocytes", labels=c("Off","On"))+
  # xlab("")+
  # ylab("Alignment Rate to hg38")+
  # ylim(0,1)+
  # theme(panel.grid.major.x = element_blank(),
		# axis.text=element_text(size=12, face="bold"),
		# axis.title=element_text(size=16) )
# dev.off()






