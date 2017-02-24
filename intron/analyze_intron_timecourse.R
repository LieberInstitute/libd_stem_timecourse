# intron expression in human stem cells
dir = '/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/'

##########################################
# load data files from intron featureCounts
load('../prep_samples/annotated_phenotype_data_stemCellTimecourse.rda')
load(paste0(dir,'intron_rpkmCounts_AZpilot_feb22.hg38.intron_n506.rda'))

#################################
# take only biological replicates
bioInd = which(pd$Class == "Naked genomes")
pd = pd[bioInd,]
pd$DAY = as.numeric(pd$DAY)
pd$DX = factor(pd$DX,levels = c('CNT','SCZ'))
pd$Donor = factor(pd$Donor)
pd$totalAssignedIntron = metrics[pd$SAMPLE_ID,'totalAssignedIntron']

intronRpkm = intronRpkm[,pd$SAMPLE_ID]

##############################
# filter intron RPKM dataframe
keepIndex = which(rowMeans(intronRpkm)>0.1)
intronRpkm = intronRpkm[keepIndex,]

#################
# make some plots
library(RColorBrewer)
library(ggplot2)
pca = prcomp(t(intronRpkm))

pdf('stem_cell_strictly_intron_timecourse_expression.pdf')
ggplot(data = data.frame(pca$x,pd))+
 geom_point(aes(x=PC1,y=PC2,fill=factor(DAY),colour = DX),pch =21)+
  scale_fill_brewer(palette = 'PiYG') + 
  scale_colour_manual(values = c('black','red'))

 ggplot(data = data.frame(pca$x,pd))+
  geom_point(aes(x=PC1,y=PC2,fill=CONDITION),pch =21)+
  scale_fill_brewer(palette = 'Set1')


ggplot(data = data.frame(pca$x,pd),aes(x=DAY,y=PC1))+
  geom_point(pch=21,aes(fill=Donor),position = position_dodge(width = .25))+
  scale_fill_brewer(palette = 'Set1') 

 ggplot(data = data.frame(pca$x,pd),aes(x=totalAssignedIntron,y=PC1))+
  geom_point(pch=21,aes(fill=Donor),position = position_dodge(width = .25))+
  scale_fill_brewer(palette = 'Set1')

 ggplot(data = data.frame(pca$x,pd),aes(x=DAY,y=PC2))+
  geom_point(pch=21,aes(fill=Donor),position = position_dodge(width = .25))+
  scale_fill_brewer(palette = 'Set1') 


pd$meanRpkm = colMeans(intronRpkm)
ggplot(data = pd,aes(x=DAY,y=meanRpkm,color = Donor))+
  geom_line()+ scale_colour_brewer(palette = 'Set1') +
  ggtitle('Mean RPKM')

ggplot(data = pd,aes(x=DAY,y=totalAssignedIntron,color = Donor))+
  geom_line()+ scale_colour_brewer(palette = 'Set1') +
  ggtitle('Intron Assignment Rate')
dev.off()


