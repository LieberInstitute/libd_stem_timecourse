library(jaffelab)
library(VariantAnnotation)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
col.pal = brewer.pal(9,"Blues")


# ##########################
# #### pd with all n = 506
# library(readxl)
# pd = read_excel("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/prep_samples/MASTER_AZ_RNA-seq_10NOV_2016.xlsx")
# pd = pd[!is.na(pd$Library),]
# pd$SampleID = paste0(pd$RNA_NO,"_",pd$Flowcell)
# ## put in order
# ord = read.csv("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/read_and_alignment_metrics_AZpilot_jan6.hg38_n506.csv", stringsAsFactors=FALSE)
# ord$RNA = ss(ord$SAMPLE_ID, "_",1)
# ord$RNA = gsub("-rerun", "", ord$RNA)
# ord$RNA = gsub("rerun", "", ord$RNA)
# ord$FC = ss(ord$SAMPLE_ID, "_",2)
# ord$id = paste0(ord$RNA,"_",ord$FC)
# pd = pd[match(ord$id, pd$SampleID),]

##########################
### pd file / sample IDs for n = 157
load("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/data/libd_stemcell_timecourse_rseGene_n157.rda")
pd = colData(rse_gene)[,1:24]

## sample indexes of 157 samples used in analyses (from 506 original samples)
id506 = read.table("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/SAMPLE_IDs.txt", stringsAsFactors=FALSE)
keepInd = which(basename(id506$V1) %in% pd$SAMPLE_ID)


##########################
### Load in snpMap
snpMap = rtracklayer::import("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Genotyping/common_missense_SNVs_hg38.bed")

##########################
# read in merged VCF file
mergedVcfFile = '/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/Genotypes/mergedVariants.vcf.gz'
vcf = readVcf(mergedVcfFile,'GRCh38.p2')
info(vcf)$RS = mcols(snpMap)$name[match(rowRanges(vcf),snpMap)]
colnames(vcf) = pd$SampleID


######################
# subset to high-depth
vcf = vcf[info(vcf)$DP > 5*ncol(vcf) & info(vcf)$DP < 80 *ncol(vcf) &
          nchar(ref(vcf)) == 1 & elementNROWS(alt(vcf)) == 1 &
          info(vcf)$VDB >0.1,]


########################################
# plot snp correlation of all samples
snps = geno(vcf)$GT
snps[snps == "."] = 0
snps[snps == "0/1"] = 1
snps[snps == "1/1"] = 2
class(snps) = "numeric"
snps = snps[,keepInd]
snpCor = cor(snps, use="pairwise.complete.obs")

rownames(snpCor) = colnames(snpCor) = paste0(pd$RNA_NO,"_",pd$Donor)
rownames(snpCor)[pd$SPECIES=="HUMAN_RAT"] = colnames(snpCor)[pd$SPECIES=="HUMAN_RAT"] = paste0(pd$RNA_NO[pd$SPECIES=="HUMAN_RAT"],"_",pd$Donor[pd$SPECIES=="HUMAN_RAT"], "_H+R")
rownames(snpCor)[pd$SPECIES=="RAT"] = colnames(snpCor)[pd$SPECIES=="RAT"] = paste0(pd$RNA_NO[pd$SPECIES=="RAT"],"_",pd$Donor[pd$SPECIES=="RAT"], "_Rat")

pdf("pheatmap_version_badoi.pdf",h=20,w=20)
pheatmap(snpCor, 
		cluster_rows=T, 
		cluster_cols=T,
		color=col.pal)
dev.off()

# pdf("pheatmap_n506.pdf",h=50,w=50)
# pheatmap(snpCor, 
		# cluster_rows=T, 
		# cluster_cols=T,
		# color=col.pal)
# dev.off()




