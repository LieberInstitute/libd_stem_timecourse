library(GenomicRanges)
library(SummarizedExperiment)
library(jaffelab)
library(ggplot2)
library(pheatmap)
theme_set(theme_bw(base_size=12) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				legend.position="none"))
library("RColorBrewer")
col.pal <- brewer.pal(9,"Blues")


### pd file / sample IDs
load("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/data/libd_stemcell_timecourse_rseGene_n157.rda")
pd = colData(rse_gene)[,2:24]
# pd = pd[which(pd$SPECIES=="HUMAN"),]
pd$ID = pd$Donor
indA = pd$ID %in% c("3","66","90")
indB = pd$ID %in% c("21","165")
pd$ID[indA] = paste0(pd$ID[indA], "_A")
pd$ID[indB] = paste0(pd$ID[indB], "_B")


### Load in variant calling files
genoFiles = paste0("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/Genotypes/", 
					pd$SAMPLE_ID, "_calledVariants.txt" )
names(genoFiles) = pd$ID
stopifnot(all(file.exists(genoFiles)))

varDat = lapply(genoFiles, function(x) {
  read.table(x, header=TRUE, fill=TRUE, sep="\t", stringsAsFactors=FALSE) })
names(varDat) = rownames(pd)  


### Load in genotype files
plink = snpStats::read.plink('AZPILOT_5_Clinical_Study_genotypes_5666161snps_imputed20140915_Straub040417')
plink_map = plink$map
plink_geno = t(as(plink$genotypes,"numeric"))

# Switch to sample ID format we're familiar with
# (based on CLIN_Sibling-Study-Fibroblasts-FBS_AUG2015_INVENTORY_14FEB2017_CUR.xlsx)
colnames(plink_geno)[grep("732",colnames(plink_geno))] = "66_A"
colnames(plink_geno)[grep("799",colnames(plink_geno))] = "3_A"
colnames(plink_geno)[grep("834",colnames(plink_geno))] = "90_A"
colnames(plink_geno)[grep("939",colnames(plink_geno))] = "21_B"
colnames(plink_geno)[grep("999164",colnames(plink_geno))] = "165_B"


## Load info on 740 SNPs used in variant calling
var_SNPs = rtracklayer::import("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Genotyping/common_missense_SNVs_hg38.bed")

## Filter genotype info to these 730 common SNPs 
var_SNPs = var_SNPs[(var_SNPs$name %in% rownames(plink_geno)), ]
plink_geno = plink_geno[var_SNPs$name,]
plink_map = plink_map[var_SNPs$name,]
plink_map$hg38_chrpos = paste0(seqnames(var_SNPs),":",start(var_SNPs))

mm = match(names(genoFiles),colnames(plink_geno))
snp = plink_geno[,mm]


pdf("hetRate.pdf") 
matchRate = lapply(seq(along=varDat), function(i) {
	x = varDat[[i]]
	s = snp[,i]
	
	x$snp_coord = paste0(x$seq, ":", x$pos)
	x$totCov = rowSums(x[,c("A","C","G","T")])
	
	m = match(x$snp_coord, plink_map$hg38_chrpos)
	x = x[!is.na(m),-10]
	s = s[m[!is.na(m)]]
	sm = plink_map[m[!is.na(m)],]
	
	## call genotypes
	#	s[x$ref == sm$allele.2] = 2-s[x$ref == sm$allele.2]
	hetRate = x$n_var/x$totCov
	hetRate[x$ref == sm$allele.2] = (x$totCov[x$ref == sm$allele.2]-x$n_var[x$ref == sm$allele.2])/x$totCov[x$ref == sm$allele.2]
	sm$highQ = x$totCov > 20 # high coverage
	sm$highQ[x$n_var < 10 & hetRate > 0.1 & hetRate < 0.9] = FALSE
	# boxplot(hetRate ~ s,ylab="HetRate",outline=FALSE,
		# main = paste(colnames(snp)[i], "- All"))
	# points(hetRate ~ jitter(s+1,amount=0.1))	
	boxplot(hetRate ~ s, subset=sm$highQ,
		main = paste0(colnames(snp)[i], " : ", pd$RNA_NO[i], ", ",pd$LibraryBatch[i],", ",pd$CONDITION[i],
					", Day ",pd$DAY[i], "\n", "Cov>20 (",pd$SPECIES[i],")" ))
	points(hetRate ~ jitter(s+1,amount=0.1),subset=sm$highQ)	

	guess = cut(hetRate, c(0,0.1,0.9,1))
	tab = table(guess,s)
	
	tab2 = table(guess[sm$highQ],s[sm$highQ])
	sm$hetRate = hetRate
	sm$genoCall = s
	sm$genoGuess = guess
	o = list(output = sm, concordAll = sum(diag(tab))/sum(tab),
		concord40 = sum(diag(tab2))/sum(tab2))
})
dev.off()
names(matchRate) = names(varDat)



ii = which.min(elementNROWS(lapply(matchRate, function(x) x$output$snp.name)))

snpNames = as.character(matchRate[[ii]]$output$snp.name)

guessSnp = sapply(matchRate, function(x) {
	x = x$output
	rownames(x) = x$snp.name
	x = x[snpNames,]
	return(as.numeric(x$genoGuess)-1)
})

## RNA-Seq Data Internal Correlation
kk = cor(guessSnp, use="pairwise.complete.obs")
rownames(kk) = paste0(pd$RNA_NO,"_",pd$Donor)
rownames(kk)[pd$SPECIES=="HUMAN_RAT"] = paste0(pd$RNA_NO[pd$SPECIES=="HUMAN_RAT"],"_",pd$Donor[pd$SPECIES=="HUMAN_RAT"], "_H+R")
colnames(kk) = paste0(pd$RNA_NO,"_",pd$Donor)
colnames(kk)[pd$SPECIES=="HUMAN_RAT"] = paste0(pd$RNA_NO[pd$SPECIES=="HUMAN_RAT"],"_",pd$Donor[pd$SPECIES=="HUMAN_RAT"], "_H+R")

pdf("pheatmap.pdf",h=20,w=20)
pheatmap(kk, 
		cluster_rows=T, 
		cluster_cols=T,
		color=col.pal)
dev.off()








































library(readxl)
library(jaffelab)
library(readr)
library(parallel)
library(ggplot2)
library(IRanges)
library(GenomicRanges)
library(pheatmap)
theme_set(theme_bw(base_size=12) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				legend.position="none"))  
library("RColorBrewer")
col.pal <- brewer.pal(9,"Blues")



## Load genotype data from the SNP-chips.
## read in genotypes of coding SNPs
setwd('/users/ssemick/BrainSwap')
load('/users/ssemick/BrainSwap/All_Platforms_Genotype_Barcode.rda')  ##"all_platforms_genotypes","Barcode_1M_map"
snpMap = Barcode_1M_map
snpMap$chromosome = gsub("23", "X", snpMap$chromosome)
snpMap$snp_coord = paste0("chr", snpMap$chromosome, ":", snpMap$position)

# split off genotypes
all_platforms_genotypes = t(all_platforms_genotypes)
chip =  ss(colnames(all_platforms_genotypes), "_",1)
chipGeno=chip
colnames(all_platforms_genotypes) = ss(colnames(all_platforms_genotypes), "_",2)
Genotypes <- all_platforms_genotypes

hg38_SNPs <- rtracklayer::import("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Genotyping/common_missense_SNVs_hg38.bed")

# SNP Chips

## Genotype Distance Matrix
## Distance matrix computed for the genotype data
## {r Genotype distance matrix}
coding_snps <- t(Genotypes[hg38_SNPs$name,])
rownames(coding_snps) <- paste0(rownames(coding_snps),"_",chipGeno)
dist_coding_snps <- as.matrix(dist(coding_snps))
#rowSums(dist_coding_snps<15)
table(rowSums(dist_coding_snps<15))


## We can test a number of different cutoffs.
## {r testing distance cutoffs}
genotype_dist_cutoff <- data.frame(sapply(1:26, function(x) {table(rowSums(dist_coding_snps<=x))[1:3]}))
genotype_dist_cutoff$Number_of_Matches <- rownames(genotype_dist_cutoff)
genotype_dist_cutoff <- tidyr::gather(genotype_dist_cutoff,key = cutoff, value = measurement, -Number_of_Matches)
genotype_dist_cutoff$cutoff <- as.numeric(gsub("X","",genotype_dist_cutoff$cutoff  ))
genotype_dist_cutoff$Number_of_Matches <- paste0(genotype_dist_cutoff$Number_of_Matches, "\n Match(s)")

c <- ggplot(genotype_dist_cutoff, aes(x=cutoff, y = measurement))
c <- c + geom_line() + facet_wrap(~Number_of_Matches, scales="free") + 
geom_vline(xintercept=15,color='red') + labs(x="Distance Threshold", y = "Number of People") 
ggsave(c,file="plots/Genotype_DistanceThreshold_Plot.pdf")
c


## Next we can compute a table of potential swaps.
## {r genodist swap table}
genoDistInd = apply(dist_coding_snps < 15, 1, which)
genoDistInd = lapply(genoDistInd, function(x) unique(names(x)))
table(elementNROWS(genoDistInd))
multi_Ind = genoDistInd[elementNROWS(genoDistInd) > 1][!duplicated(genoDistInd[elementNROWS(genoDistInd) > 1])]

genoDistInd2 = genoDistInd[elementNROWS(genoDistInd) == 2]
genoDistInd3 = genoDistInd[elementNROWS(genoDistInd) == 3]

matched2 = unlist(genoDistInd2[!duplicated(genoDistInd2)])
wrong_matched2<-matched2[(!gsub("\\_.*","",matched2) == gsub("\\_.*","",names(matched2)))]

genoDistInd3[!duplicated(genoDistInd3)]
#matched3 = unlist(genoDistInd3[!duplicated(genoDistInd3)])
#wrong_matched3 <- matched3[(!gsub("\\_.*","",matched3) == gsub("\\_.*","",names(matched3)))]
wrong3 <-t(sapply(genoDistInd3[!duplicated(genoDistInd3)],c))
colnames(wrong3) <-c("BrNum1","BrNum2","BrNum3")

Geno_misMatch = data.frame(BrNum1 = names(wrong_matched2), BrNum2 = wrong_matched2, BrNum3=NA)
Geno_misMatch <- rbind(Geno_misMatch,wrong3)
Geno_misMatch[c('BrNum1','Platform1')] = c(ss(as.character(Geno_misMatch$BrNum1),"_",1), ss(as.character(Geno_misMatch$BrNum1),"_",2) )
Geno_misMatch[c('BrNum2','Platform2')] = c(ss(as.character(Geno_misMatch$BrNum2),"_",1), ss(as.character(Geno_misMatch$BrNum2),"_",2) )
Geno_misMatch[c('BrNum3','Platform3')] = c(ss(as.character(Geno_misMatch$BrNum3),"_",1), ss(as.character(Geno_misMatch$BrNum3),"_",2) )

Geno_misMatch$Platform1[is.na(Geno_misMatch$Platform3)] <- substr(Geno_misMatch$Platform1[is.na(Geno_misMatch$Platform3)],1,nchar(Geno_misMatch$Platform1[is.na(Geno_misMatch$Platform3)])-1)

Geno_misMatch[c('BrNum1','BrNum2','Platform1','Platform2')] <- t(apply(Geno_misMatch[,c('BrNum1','BrNum2','Platform1','Platform2')],1, function(x) {id=order(x[1:2])
c(x[1:2][id],x[3:4][id])}) )

knitr::kable(Geno_misMatch,row.names=FALSE)
```

## Genotype Correlation Matrix
## Correlation matrix computed for the genotype data
## {r Genotype correlation matrix}
tmp <- all_platforms_genotypes
colnames(tmp) <- paste0(colnames(tmp),"_",chipGeno)
cor_coding_snps <- cor(tmp[,order(gsub("_.*","",colnames(tmp)))], use="pairwise.complete.obs")
#pdf('plots/Genotype_Cor_pheatmap.pdf')
#pheatmap(cor_coding_snps, cluster_rows=FALSE, cluster_cols=FALSE)
#dev.off()
```

## We can test a number of different cutoffs.
## {r testing cor cutoffs}
genotype_cor_cutoff <- data.frame(sapply(seq(from = 0.45, to = .99, by = 0.001), function(x) {table(rowSums(cor_coding_snps>=x))[1:3]}))
colnames(genotype_cor_cutoff) <- seq(from = 0.45, to = .99, by = 0.001)

genotype_cor_cutoff$Number_of_Matches <- 1:3
genotype_cor_cutoff <- tidyr::gather(genotype_cor_cutoff,key = cutoff, value = measurement, -Number_of_Matches)
genotype_cor_cutoff$cutoff <- as.numeric(genotype_cor_cutoff$cutoff)
genotype_cor_cutoff$Number_of_Matches <- paste0(genotype_cor_cutoff$Number_of_Matches, "\n Match(s)")

d <- ggplot(genotype_cor_cutoff, aes(x=cutoff, y = measurement))
d <- d +  geom_line() + facet_wrap(~Number_of_Matches, scales="free") +
geom_vline(xintercept=0.6,color='red') + 
labs(x="Correlation Threshold", y = "Number of People") 
ggsave(d,file="plots/Genotype_PearsonCorThreshold_Plot.pdf")
d
```

## {r genocor swap table}
genoCorInd = apply(cor_coding_snps > 0.6, 1, which)
genoCorInd = lapply(genoCorInd, function(x) unique(names(x)))
table(elementNROWS(genoCorInd))
multi_Ind = genoCorInd[elementNROWS(genoCorInd) > 1][!duplicated(genoCorInd[elementNROWS(genoCorInd) > 1])]

genoCorInd2 = genoCorInd[elementNROWS(genoCorInd) == 2]
genoCorInd3 = genoCorInd[elementNROWS(genoCorInd) == 3]

matched2 = unlist(genoCorInd2[!duplicated(genoCorInd2)])
wrong_matched2<-matched2[(!gsub("\\_.*","",matched2) == gsub("\\_.*","",names(matched2)))]


genoCorInd3[!duplicated(genoCorInd3)]
#matched3 = unlist(genoCorInd3[!duplicated(genoCorInd3)])
wrong3 <-t(sapply(genoCorInd3[!duplicated(genoCorInd3)],c))
colnames(wrong3) <-c("BrNum1","BrNum2","BrNum3")

Geno_misMatch = data.frame(BrNum1 = names(wrong_matched2), BrNum2 = wrong_matched2, BrNum3=NA)
Geno_misMatch <- rbind(Geno_misMatch,wrong3)
Geno_misMatch[c('BrNum1','Platform1')] = c(ss(as.character(Geno_misMatch$BrNum1),"_",1), ss(as.character(Geno_misMatch$BrNum1),"_",2) )
Geno_misMatch[c('BrNum2','Platform2')] = c(ss(as.character(Geno_misMatch$BrNum2),"_",1), ss(as.character(Geno_misMatch$BrNum2),"_",2) )
Geno_misMatch[c('BrNum3','Platform3')] = c(ss(as.character(Geno_misMatch$BrNum3),"_",1), ss(as.character(Geno_misMatch$BrNum3),"_",2) )

Geno_misMatch$Platform1[is.na(Geno_misMatch$Platform3)] <- substr(Geno_misMatch$Platform1[is.na(Geno_misMatch$Platform3)],1,nchar(Geno_misMatch$Platform1[is.na(Geno_misMatch$Platform3)])-1)

knitr::kable(Geno_misMatch,row.names=FALSE)
```

## {r plotting genotype correlation swap}
genoSwap_involved <- unique(c(Geno_misMatch$BrNum1,Geno_misMatch$BrNum2,Geno_misMatch$BrNum3))
genoSwap_involved[!is.na(genoSwap_involved)]
Geno_corSwap = cor_coding_snps[gsub("_.*","",rownames(cor_coding_snps)) %in% genoSwap_involved, gsub("_.*","",rownames(cor_coding_snps)) %in% genoSwap_involved ]

Geno_corSwap_anno <- data.frame(Chip = gsub(".*_","", rownames(Geno_corSwap) ), row.names = rownames(Geno_corSwap) )

##### Heatmap #####
pdf('/users/eburke/Genotype_Cor_pheatmap_corSWAP.pdf', height=11, width=8.5)
pheatmap(Geno_corSwap, 
		cluster_rows=T, 
		cluster_cols=T, 
		labels_row = gsub("_.*","",rownames(Geno_corSwap)), 
		labels_col = gsub("_.*","",colnames(Geno_corSwap)),
		annotation_row = Geno_corSwap_anno,
		color=col.pal)
pheatmap(Geno_corSwap, 
		cluster_rows=T, 
		cluster_cols=T, 
		labels_row = gsub("_.*","",rownames(Geno_corSwap)), 
		labels_col = gsub("_.*","",colnames(Geno_corSwap)),
		annotation_row = Geno_corSwap_anno,
		color = col.pal)
dev.off()


















# DLPFC RiboZero 

## Loading Phenotype Data
Load phenotype data from dlpfc RiboZero Cohort.
```{r load DLPFC RiboZero}
pd <- read.delim('/users/ajaffe/Lieber/Projects/RNAseq/Consortium/PhaseII_DLPFC/DLPFC_PhaseII_sample-selection_04_20_2015.txt', header=TRUE)

colnames(pd)[colnames(pd)=="BRNum"] <-"BrNum"
pd = pd[!is.na(pd$BrNum),]
pd = as.data.frame(pd)
pd$BrNum = paste0("Br", as.numeric(pd$BrNum))
pd$RNum = paste0("R", as.numeric(pd$RNum))

pd$hasGeno = pd$BrNum %in% colnames(Genotypes)
pd = pd[pd$hasGeno,] ## filter for now

pd <- pd[!duplicated(pd$RNum),]
```

## Filtering
Filter down to sequenced subjects with genotype data.
```{r filter DLPFC RiboZero}
########################################
## filter subjects w/o genotype data ###
########################################

## read in genotypes 
fn = list.files("/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/Genotypes",full.names=TRUE)
names(fn) = ss( list.files("/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/Genotypes"), "_")
flow = ss( list.files("/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/Genotypes"), "_", 2)
names(fn) = pd$BrNum[match(names(fn), pd$RNum)]
flow = flow[!is.na(names(fn))]
fn = fn[!is.na(names(fn))]

## read in
#genoDat = mclapply(fn, read.delim, as.is=TRUE,mc.cores=8)
genoDat = mclapply(fn, read.delim, as.is=TRUE,mc.cores=1)
```

## hetRate Calculation
Calculate `hetRate` and make plots for the samples.
```{r hetRate DLPFC RiboZero}
Genotypes <- Genotypes[hg38_SNPs$name,]
snpMap <- snpMap[hg38_SNPs$name,]
snpMap$hg38_chrpos <-paste0(seqnames(hg38_SNPs),":",start(hg38_SNPs))
 
mm = match(names(fn),colnames(Genotypes))
snp = Genotypes[,mm]
chip = chipGeno[mm]

pdf("plots/hetRate_check_DLPFC_RiboZero.pdf",h=4,w=8)
mypar(1,2)
matchRate = lapply(seq(along=genoDat), function(i) {
#	cat(".")
	x = genoDat[[i]]
	s = snp[,i]
	
	x$snp_coord = paste0(x$seq, ":", x$pos)
	x$totCov = rowSums(x[,c("A","C","G","T")])
	
	m = match(x$snp_coord, snpMap$hg38_chrpos)
	x = x[!is.na(m),-10]
	s = s[m[!is.na(m)]]
	sm = snpMap[m[!is.na(m)],]
	
	## call genotypes
	#	s[x$ref == sm$allele.2] = 2-s[x$ref == sm$allele.2]
	hetRate = x$n_var/x$totCov
	hetRate[x$ref == sm$allele.2] = (x$totCov[x$ref == sm$allele.2]-x$n_var[x$ref == sm$allele.2])/x$totCov[x$ref == sm$allele.2]
	sm$highQ = x$totCov > 20 # high coverage
	sm$highQ[x$n_var < 10 & hetRate > 0.1 & hetRate < 0.9] = FALSE
	boxplot(hetRate ~ s,ylab="HetRate",outline=FALSE,
		main = paste(colnames(snp)[i], "- All"))
	points(hetRate ~ jitter(s+1,amount=0.1))	
	boxplot(hetRate ~ s, subset=sm$highQ,
		main = paste(colnames(snp)[i], "- Cov>20"))
	points(hetRate ~ jitter(s+1,amount=0.1),subset=sm$highQ)	

	guess = cut(hetRate, c(0,0.1,0.9,1))
	tab = table(guess,s)
	
	tab2 = table(guess[sm$highQ],s[sm$highQ])
	sm$hetRate = hetRate
	sm$genoCall = s
	sm$genoGuess = guess
	o = list(output = sm, concordAll = sum(diag(tab))/sum(tab),
		concord40 = sum(diag(tab2))/sum(tab2))
})
dev.off()
names(matchRate) = names(genoDat)
save(matchRate, file="rdas/matchOutput_DLPFC_RiboZero.rda")
```

Next we checked the `matchRate` with and without a quality filter.
```{r calculate matchRate DLPFC RiboZero}
pd <- pd[match(names(matchRate),pd$BrNum),]
pd$matchRate = sapply(matchRate, function(x) x$concord40)
pd$matchRate_All = sapply(matchRate, function(x) x$concordAll)
cor.test(pd$RIN,pd$matchRate)
```

Matchrate Histogram
```{r histo with the quality filter DLPFC RiboZero}
dat2 <- tidyr::gather(pd[,c('matchRate','matchRate_All')])

matchRate_histo <- ggplot(data = dat2, aes(x=value))
matchRate_histo <- matchRate_histo + geom_histogram( binwidth=0.01, 
								linetype=1, 
								alpha=.25, 
								color="black",
								fill="red") +
								labs(x="Match Rate") +
								scale_y_continuous() + 
								scale_x_continuous() + facet_wrap(~key)
ggsave('plots/matchRate_histo_DLPFC_RiboZero.pdf',height=8.5,width=11)
matchRate_histo
```

These are the `BrNums` that have a `matchRate` below 0.8.
```{r poor match rate DLPFC RiboZero}
unique(pd$BrNum[which(pd$matchRate < 0.8)])
pd[which(pd$matchRate < 0.8),]
```

## Consensus Set of SNPs
```{r consensus set of SNPs: guessing DLPFC RiboZero}
ii = which.min(elementNROWS(lapply(matchRate, function(x) x$output$snp.name)))

snpNames = as.character(matchRate[[ii]]$output$snp.name)

guessSnp = sapply(matchRate, function(x) {
	x = x$output
	rownames(x) = x$snp.name
	x = x[snpNames,]
	return(as.numeric(x$genoGuess)-1)
})
```

## RNA-Seq Data Internal Correlation
```{r correlation of guessed SNP DLPFC RiboZero}
kk = cor(guessSnp, use="pairwise.complete.obs")
corInd = apply(kk > 0.6, 1, which)
corInd = lapply(corInd, function(x) unique(names(x)))
table(elementNROWS(corInd))

pheatmap(kk, 
		cluster_rows=T, 
		cluster_cols=T,
		color=col.pal)
```
Each person matches uniquely thus there is no evidence of sample duplication within the RNA-Seq data.

```{r key objects DLPFC RiboZero}
matchRate_DLPFC_RiboZero <- matchRate
flow_DLPFC_RiboZero <- flow
pd_DLPFC_RiboZero <- pd
```

# Hippo RiboZero PhaseII 

## Loading Phenotype Data
Load phenotype data from HippoRiboZero Cohort.
```{r loading Hippo RiboZero}
# phenotype data
# original path: /users/ajaffe/Lieber/Projects/RNAseq/Consortium/PhaseII_HIPPO/PhaseII_Hippo_LIBD_finalSet.txt
pd <- read.delim('/users/ajaffe/Lieber/Projects/RNAseq/Consortium/PhaseII_HIPPO/PhaseII_Hippo_LIBD_finalSet.txt', header = TRUE)
colnames(pd)[colnames(pd)=="BRNum"] <-"BrNum"
pd = pd[!is.na(pd$BrNum),]
pd = as.data.frame(pd)
pd$BrNum = paste0("Br", as.numeric(pd$BrNum))
pd$RnumHIPPOPhaseII = paste0("R", as.numeric(pd$RnumHIPPOPhaseII))

pd$hasGeno = pd$BrNum %in% colnames(Genotypes)
pd = pd[pd$hasGeno,] ## filter for now

pd <- pd[!duplicated(pd$RnumHIPPOPhaseII),]
```

## Filtering
Filter down to sequenced subjects with genotype data.
```{r filter HIPPO RiboZero}
########################################
## filter subjects w/o genotype data ###
########################################

## read in genotypes 
fn = list.files("/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/Genotypes",full.names=TRUE)
names(fn) = ss( list.files("/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/Genotypes"), "_")
flow = ss( list.files("/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/Genotypes"), "_", 2)
names(fn) = pd$BrNum[match(names(fn), pd$RnumHIPPOPhaseII)]
flow = flow[!is.na(names(fn))]
fn = fn[!is.na(names(fn))]

## read in
#genoDat = mclapply(fn, read.delim, as.is=TRUE,mc.cores=8)
genoDat = mclapply(fn, read.delim, as.is=TRUE,mc.cores=1)
```

## hetRate Calculation
Calculate `hetRate` and make plots for the samples.
```{r hetRate Hippo RiboZero}
Genotypes <- Genotypes[hg38_SNPs$name,]
snpMap <- snpMap[hg38_SNPs$name,]
snpMap$hg38_chrpos <-paste0(seqnames(hg38_SNPs),":",start(hg38_SNPs))

mm = match(names(fn),colnames(Genotypes))
snp = Genotypes[,mm]
chip = chipGeno[mm]

pdf("plots/hetRate_check_HIPPO_RiboZero.pdf",h=4,w=8)
mypar(1,2)
matchRate = lapply(seq(along=genoDat), function(i) {
#	cat(".")
	x = genoDat[[i]]
	s = snp[,i]
	
	x$snp_coord = paste0(x$seq, ":", x$pos)
	x$totCov = rowSums(x[,c("A","C","G","T")])
	
	m = match(x$snp_coord, snpMap$hg38_chrpos)
	x = x[!is.na(m),-10]
	s = s[m[!is.na(m)]]
	sm = snpMap[m[!is.na(m)],]
	
	## call genotypes
	#	s[x$ref == sm$allele.2] = 2-s[x$ref == sm$allele.2]
	hetRate = x$n_var/x$totCov
	hetRate[x$ref == sm$allele.2] = (x$totCov[x$ref == sm$allele.2]-x$n_var[x$ref == sm$allele.2])/x$totCov[x$ref == sm$allele.2]
	sm$highQ = x$totCov > 20 # high coverage
	sm$highQ[x$n_var < 10 & hetRate > 0.1 & hetRate < 0.9] = FALSE
	boxplot(hetRate ~ s,ylab="HetRate",outline=FALSE,
		main = paste(colnames(snp)[i], "- All"))
	points(hetRate ~ jitter(s+1,amount=0.1))	
	boxplot(hetRate ~ s, subset=sm$highQ,
		main = paste(colnames(snp)[i], "- Cov>20"))
	points(hetRate ~ jitter(s+1,amount=0.1),subset=sm$highQ)	

	guess = cut(hetRate, c(0,0.1,0.9,1))
	tab = table(guess,s)
	
	tab2 = table(guess[sm$highQ],s[sm$highQ])
	sm$hetRate = hetRate
	sm$genoCall = s
	sm$genoGuess = guess
	o = list(output = sm, concordAll = sum(diag(tab))/sum(tab),
		concord40 = sum(diag(tab2))/sum(tab2))
})
dev.off()
names(matchRate) = names(genoDat)
save(matchRate, file="rdas/matchOutput_HIPPO_RiboZero.rda")
```

Next we checked the `matchRate` with and without a quality filter.
```{r calculate matchRate Hippo RiboZero}
pd <- pd[match(names(matchRate),pd$BrNum),]
pd$matchRate = sapply(matchRate, function(x) x$concord40)
pd$matchRate_All = sapply(matchRate, function(x) x$concordAll)
cor.test(pd$RIN,pd$matchRate)
```

Matchrate Histogram
```{r histo with the quality filter Hippo RiboZero}
dat2 <- tidyr::gather(pd[,c('matchRate','matchRate_All')])

matchRate_histo <- ggplot(data = dat2, aes(x=value))
matchRate_histo <- matchRate_histo + geom_histogram( binwidth=0.01, 
								linetype=1, 
								alpha=.25, 
								color="black",
								fill="red") +
								labs(x="Match Rate") +
								scale_y_continuous() + 
								scale_x_continuous() + facet_wrap(~key)
ggsave('plots/matchRate_histo_HIPPO_RiboZero.pdf',height=8.5,width=11)
matchRate_histo
```

These are the `BrNums` that have a `matchRate` below 0.8.
```{r poor match rate Hippo RiboZero}
unique(pd$BrNum[which(pd$matchRate < 0.8)])
pd[which(pd$matchRate < 0.8),]
```

## Consensus Set of SNPs
```{r consensus set of SNPs: guessing Hippo RiboZero}
ii = which.min(elementNROWS(lapply(matchRate, function(x) x$output$snp.name)))

snpNames = as.character(matchRate[[ii]]$output$snp.name)

guessSnp = sapply(matchRate, function(x) {
	x = x$output
	rownames(x) = x$snp.name
	x = x[snpNames,]
	return(as.numeric(x$genoGuess)-1)
})
```

## RNA-Seq Data Internal Correlation
```{r correlation of guessed SNP Hippo RiboZero}
kk = cor(guessSnp, use="pairwise.complete.obs")
corInd = apply(kk > 0.6, 1, which)
corInd = lapply(corInd, function(x) unique(names(x)))
table(elementNROWS(corInd))

pheatmap(kk, 
		cluster_rows=T, 
		cluster_cols=T,
		color=col.pal)
```
Each person matches uniquely thus there is no evidence of sample duplication within the RNA-Seq data.

```{r key objects Hippo RiboZero}
matchRate_HIPPO_RiboZero <- matchRate
flow_HIPPO_RiboZero <- flow
pd_HIPPO_RiboZero <- pd
```

# Dentate Gyrus (Paired End)

## Loading Phenotype Data
Load phenotype data from Dentate Gyrus PE Cohort.
```{r loading phenotype data DGPE}
# phenotype data
# original path: /users/ajaffe/Lieber/Projects/RNAseq/Consortium/PhaseII_HIPPO/PhaseII_Hippo_LIBD_finalSet.txt
pd <- read.csv('/dcl01/lieber/ajaffe/lab/dg_hippo/preprocessed_data/pd.csv', header = TRUE)
colnames(pd)[colnames(pd)=="BRNum"] <-"BrNum"
pd = pd[!is.na(pd$BrNum),]
pd = as.data.frame(pd)

pd$hasGeno = pd$BrNum %in% colnames(Genotypes)
pd = pd[pd$hasGeno,] ## filter for now

pd <-pd[pd$Prep=="paired_end",] #filter to paired end only

```

## Filtering
Filter down to sequenced subjects with genotype data.
```{r filter DGPE PhaseIII}
########################################
## filter subjects w/o genotype data ###
########################################

## read in genotypes 
fn = list.files("/dcl01/lieber/ajaffe/lab/dg_hippo/preprocessed_data/paired_end_n292/Genotypes",full.names=TRUE)
names(fn) = ss( list.files("/dcl01/lieber/ajaffe/lab/dg_hippo/preprocessed_data/paired_end_n292/Genotypes"), "_")
flow = ss( list.files("/dcl01/lieber/ajaffe/lab/dg_hippo/preprocessed_data/paired_end_n292/Genotypes"), "_", 2)
names(fn) = pd$BrNum[match(names(fn), pd$RNum)]
flow = flow[!is.na(names(fn))]
flow = gsub("calledVariants.txt", "", flow)
fn = fn[!is.na(names(fn))]

## read in
#genoDat = mclapply(fn, read.delim, as.is=TRUE,mc.cores=8)
genoDat = mclapply(fn, read.delim, as.is=TRUE,mc.cores=1)
```

## hetRate Calculation
Calculate `hetRate` and make plots for the samples.
```{r hetRate calculations DGPE PhaseIII}
Genotypes <- Genotypes[hg38_SNPs$name,]
snpMap <- snpMap[hg38_SNPs$name,]
snpMap$hg38_chrpos <-paste0(seqnames(hg38_SNPs),":",start(hg38_SNPs))

mm = match(names(fn),colnames(Genotypes))
snp = Genotypes[,mm]
chip = chipGeno[mm]

pdf("plots/hetRate_check_DentateGyrus_PE.pdf",h=4,w=8)
mypar(1,2)
matchRate = lapply(seq(along=genoDat), function(i) {
#	cat(".")
	x = genoDat[[i]]
	s = snp[,i]
	
	x$snp_coord = paste0(x$seq, ":", x$pos)
	x$totCov = rowSums(x[,c("A","C","G","T")])
	
	m = match(x$snp_coord, snpMap$hg38_chrpos)
	x = x[!is.na(m),-10]
	s = s[m[!is.na(m)]]
	sm = snpMap[m[!is.na(m)],]
	
	## call genotypes
	#	s[x$ref == sm$allele.2] = 2-s[x$ref == sm$allele.2]
	hetRate = x$n_var/x$totCov
	hetRate[x$ref == sm$allele.2] = (x$totCov[x$ref == sm$allele.2]-x$n_var[x$ref == sm$allele.2])/x$totCov[x$ref == sm$allele.2]
	sm$highQ = x$totCov > 20 # high coverage
	sm$highQ[x$n_var < 10 & hetRate > 0.1 & hetRate < 0.9] = FALSE
	boxplot(hetRate ~ s,ylab="HetRate",outline=FALSE,
		main = paste(colnames(snp)[i], "- All"))
	points(hetRate ~ jitter(s+1,amount=0.1))	
	boxplot(hetRate ~ s, subset=sm$highQ,
		main = paste(colnames(snp)[i], "- Cov>20"))
	points(hetRate ~ jitter(s+1,amount=0.1),subset=sm$highQ)	

	guess = cut(hetRate, c(0,0.1,0.9,1))
	tab = table(guess,s)
	
	tab2 = table(guess[sm$highQ],s[sm$highQ])
	sm$hetRate = hetRate
	sm$genoCall = s
	sm$genoGuess = guess
	o = list(output = sm, concordAll = sum(diag(tab))/sum(tab),
		concord40 = sum(diag(tab2))/sum(tab2))
})
dev.off()
names(matchRate) = names(genoDat)
save(matchRate, file="rdas/matchOutput_DentateGyrus_PE.rda")
```

Next we checked the `matchRate` with and without a quality filter.
```{r calculate matchRate DGPE PhaseIII}
pd <- pd[match(names(matchRate),pd$BrNum),]
pd$matchRate = sapply(matchRate, function(x) x$concord40)
pd$matchRate_All = sapply(matchRate, function(x) x$concordAll)
cor.test(pd$RIN,pd$matchRate)
```

Matchrate Histogram
```{r histo with the quality filter DGPE PhaseIII}
dat2 <- tidyr::gather(pd[,c('matchRate','matchRate_All')])

matchRate_histo <- ggplot(data = dat2, aes(x=value))
matchRate_histo <- matchRate_histo + geom_histogram( binwidth=0.01, 
								linetype=1, 
								alpha=.25, 
								color="black",
								fill="red") +
								labs(x="Match Rate") +
								scale_y_continuous() + 
								scale_x_continuous() + facet_wrap(~key)
ggsave('plots/matchRate_histo_DGEPE.pdf',height=8.5,width=11)
matchRate_histo
```

These are the `BrNums` that have a `matchRate` below 0.8.
```{r poor match rate DGPE PhaseIII}
unique(pd$BrNum[which(pd$matchRate < 0.8)])
pd[which(pd$matchRate < 0.8),]
```

## Consensus Set of SNPs
```{r consensus set of SNPs: guessing DGPE PhaseIII}
ii = which.min(elementNROWS(lapply(matchRate, function(x) x$output$snp.name)))

snpNames = as.character(matchRate[[ii]]$output$snp.name)

guessSnp = sapply(matchRate, function(x) {
	x = x$output
	rownames(x) = x$snp.name
	x = x[snpNames,]
	return(as.numeric(x$genoGuess)-1)
})
```

## RNA-Seq Data Internal Correlation
```{r correlation of guessed SNP DGPE PhaseIII}
kk = cor(guessSnp, use="pairwise.complete.obs")
corInd = apply(kk > 0.6, 1, which)
corInd = lapply(corInd, function(x) unique(names(x)))
table(elementNROWS(corInd))

pheatmap(kk, 
		cluster_rows=T, 
		cluster_cols=T,
		color=col.pal)
```
Each person matches uniquely thus there is no evidence of sample duplication within the RNA-Seq data.

```{r key objects DGPE PhaseIII}
matchRate_DGPE <- matchRate
flow_DGPE <- flow
pd_DGPE <- pd
```

```{r extracting common SNPs across all platforms}
#HIPPO RiboZero 
ii = which.min(elementNROWS(lapply(matchRate_HIPPO_RiboZero, function(x) x$output$snp.name)))
snpNames_HIPPO_RiboZero = as.character(matchRate_HIPPO_RiboZero[[ii]]$output$snp.name)
#DLPFC RiboZero 
ii = which.min(elementNROWS(lapply(matchRate_DLPFC_RiboZero, function(x) x$output$snp.name)))
snpNames_DLPFC_RiboZero = as.character(matchRate_DLPFC_RiboZero[[ii]]$output$snp.name)
#Dentate Gyrus Paired End 
ii = which.min(elementNROWS(lapply(matchRate_DGPE, function(x) x$output$snp.name)))
snpNames_DGPE = as.character(matchRate_DGPE[[ii]]$output$snp.name)

# Commmon Snp Names
snpNames <- Reduce(intersect, list(snpNames_HIPPO_RiboZero,
                                   snpNames_DLPFC_RiboZero,
                                   snpNames_DGPE))

guessSnp = function(matchRate) { sapply(matchRate, function(x) {
	x = x$output
	rownames(x) = x$snp.name
	x = x[snpNames,]
	return(as.numeric(x$genoGuess)-1)
}) }

snp_HIPPO_RiboZero <- guessSnp(matchRate = matchRate_HIPPO_RiboZero)
colnames(snp_HIPPO_RiboZero) <-  paste0(colnames(snp_HIPPO_RiboZero),"_HippoRiboZero")
colnames(snp_HIPPO_RiboZero)<-paste0(colnames(snp_HIPPO_RiboZero),":",flow_HIPPO_RiboZero)

snp_DLPFC_RiboZero <- guessSnp(matchRate_DLPFC_RiboZero)
colnames(snp_DLPFC_RiboZero) <-  paste0(colnames(snp_DLPFC_RiboZero),"_DlpfcRiboZero")
colnames(snp_DLPFC_RiboZero)<-paste0(colnames(snp_DLPFC_RiboZero),":",flow_DLPFC_RiboZero)

snp_DGEPE <- guessSnp(matchRate_DGPE)
colnames(snp_DGEPE) <-  paste0(colnames(snp_DGEPE),"_DentateGyrus")
colnames(snp_DGEPE)<-paste0(colnames(snp_DGEPE),":",flow_DGPE)
```

## Matching on All People Between Genotypes and RNA-seq
```{r consensus set of SNPs All People: calling}
colnames(Genotypes) <- paste0(colnames(Genotypes),"_", chipGeno)
```

```{r merging RNAseq with genotype data}
RNAseq_Geno <- cbind(Genotypes[snpNames,],snp_DLPFC_RiboZero,snp_HIPPO_RiboZero, snp_DGEPE)
cor_all_methods <- cor(RNAseq_Geno[,order(gsub("_.*","",colnames(RNAseq_Geno)))], use="pairwise.complete.obs")
```

We can test a number of different cutoffs.
```{r testing cor cutoffs for all}
matchMax <- length(table(rowSums(cor_all_methods>=.6)))
test_seq <- seq(from = 0.35, to = .95, by = 0.05)
all_cor_cutoff <- data.frame(sapply(test_seq, function(x) {table(rowSums(cor_all_methods>=x))[1:matchMax]}))
colnames(all_cor_cutoff) <- test_seq

all_cor_cutoff$Number_of_Matches <- 1:matchMax
all_cor_cutoff <- tidyr::gather(all_cor_cutoff,key = cutoff, value = measurement, -Number_of_Matches)
all_cor_cutoff$cutoff <- as.numeric(all_cor_cutoff$cutoff)
all_cor_cutoff$Number_of_Matches <- paste0(all_cor_cutoff$Number_of_Matches, "\n Match(s)")

d <- ggplot(all_cor_cutoff, aes(x=cutoff, y = measurement))
d <- d +  geom_line() + facet_wrap(~Number_of_Matches, scales="free") +
geom_vline(xintercept=0.6,color='red') + 
labs(x="Correlation Threshold", y = "Number of People") 
ggsave(d,file="plots/All_RNASeq_Genotype_PearsonCorThreshold_Plot.pdf", height=8.5,width=11)
d
```

```{r ALL swap table}
allCorInd = apply(cor_all_methods > 0.6, 1, which)
allCorInd = lapply(allCorInd, function(x) unique(names(x)))
table(elementNROWS(allCorInd))
multi_Ind = allCorInd[elementNROWS(allCorInd) > 1][!duplicated(allCorInd[elementNROWS(allCorInd) > 1])]

indx <- lengths(multi_Ind) 
 res <- as.data.frame(do.call(rbind,lapply(multi_Ind, `length<-`,
                          max(indx))))
colnames(res) <- paste0("BrNum",1:ncol(res) )
res<-sapply(res,as.character)
res_wrong <- res[ which(apply((gsub("_.*","",res)), 1, function(x) { x<-x[!is.na(x)] 
sum (duplicated(x) | duplicated(x, fromLast=TRUE))!=length(x) }) ),] #pick out rows where there is at least ONE disagreement about genotype and brnum
res_wrong <- res_wrong[,-which(colSums(is.na(res_wrong))==nrow(res_wrong))] #drop columns with all NA

#matched2 = unlist(allCorInd2[!duplicated(allCorInd2)])

Geno_misMatch = as.data.frame(res_wrong[order(rowSums(is.na(res_wrong)),decreasing=F ),])
brnum <- sapply(Geno_misMatch, function(x) c(ss(as.character(x),"_",1) ) )
platform <- sapply(Geno_misMatch, function(x) c(ss(as.character(x),"_",2) ) )
colnames(platform) <- gsub("BrNum","Platform",colnames(platform))
Geno_misMatch <- cbind(brnum,platform)

knitr::kable(Geno_misMatch,row.names=FALSE)
```

```{r plotting genocor ALL correlation swap}
genoSwap_involved <- unique(Geno_misMatch[,grepl("BrNum",colnames(Geno_misMatch))] )
genoSwap_involved<-genoSwap_involved[!is.na(genoSwap_involved)]
All_corSwap = cor_all_methods[gsub("_.*","",rownames(cor_all_methods)) %in% genoSwap_involved, gsub("_.*","",rownames(cor_all_methods)) %in% genoSwap_involved ]

All_corSwap_anno <- data.frame(Chip = gsub(".*_","", rownames(All_corSwap) ), row.names = rownames(All_corSwap) )
All_corSwap_anno$Chip <- gsub(":.*","",All_corSwap_anno$Chip)

pdf('plots/All_RNAseq_All_Cor_pheatmap_corSWAP.pdf', height=11, width=8.5,onefile = FALSE)
pheatmap(All_corSwap, 
		cluster_rows=T, 
		cluster_cols=T, 
		labels_row = gsub("_.*","",rownames(All_corSwap)), 
		labels_col = gsub("_.*","",colnames(All_corSwap)),
		annotation_row = All_corSwap_anno,
		color=col.pal)
dev.off()
pheatmap(All_corSwap, 
		cluster_rows=T, 
		cluster_cols=T, 
		labels_row = gsub("_.*","",rownames(All_corSwap)), 
		labels_col = gsub("_.*","",colnames(All_corSwap)),
		annotation_row = All_corSwap_anno,
		color=col.pal)
# fn[duplicated(paste0(names(fn),flow))|duplicated(paste0(names(fn),flow), fromLast=TRUE)]

```
