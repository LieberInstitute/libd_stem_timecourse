###

library(jaffelab)
library(readxl)
library(lubridate)

## read in biological phenotype info
pheno = read_excel("MASTER_AZ_RNA-seq_10NOV_2016.xlsx")

## rename columns, some have changed
names(pheno)[c(9:10,12,16:28,31,43)] = c("LibConstruct",
	"DateToAJ", "JC_Comment", "NumReads", "SeqNotes",
	"NumMapped", "mapRate", "NumProperMap", "properRate",
	"NumUnmapped", "unmapRate",
	"riboMapped", "riboRate", "mitoMapped", "mitoRate",
	"alignNote", "LINE","Fluidigm")

## drop old tophat alignment info
pheno = pheno[,c(4:15, 29:37, 39:ncol(pheno))]

## make RNA Number and flowcell the unique ID
pheno$SampleID = paste0(pheno$RNA_NO, "_", pheno$Flowcell)

#####################
## read in alignment info to select samples
pd = read.csv("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Projects/AZpilot_stemcell/read_and_alignment_metrics_AZpilot_jan6.hg38_n506.csv",
	as.is=TRUE,row.names=1)
pd$SampleID = paste0(ss(pd$SAMPLE_ID, "_", 1), "_", 
	ss(pd$SAMPLE_ID, "_", 2))

## match up
pd = merge(pheno, pd, all=FALSE, by="SampleID")
rownames(pd) = pd$SampleID

###################
## filter the phenotype data to timecourse 
##	and technical samples 
timecourseIndex = which(pd$Class %in% c("Naked genomes",
		"Internal Control") &
	pd$SPECIES %in% c("HUMAN", "HUMAN_RAT"))
pd = pd[timecourseIndex,]

## which have fluidigm
dat = read_excel("../fluidigm/FLUIDIGM_IPS_01252016.xlsx",	sheet=1)
colnames(dat)[1] = "Sample"

table(pd$LINE %in% dat$Sample)
table(pd$LINE[!pd$LINE %in% dat$Sample])

############################
## clean up some columns

## batch
pd$LibConstruct = ymd(pd$LibConstruct)

pd$Library[pd$Library == "BATCH10"] = "Batch10"
colnames(pd)[6] = "LibraryBatch"
pd$LibraryBatch = factor(pd$LibraryBatch, 
	levels = paste0("Batch", 1:10))
pd$LibraryBatch  = droplevels(pd$LibraryBatch )

## drop unused columns
pd = pd[,-(19:23)]

save(pd, file="annotated_phenotype_data_stemCellTimecourse.rda"
## 

## keep low map rate for now...
boxplot(overallMapRate ~ CONDITION, data=pd)
boxplot(totalAssignedGene ~ CONDITION, data=pd)
plot(totalAssignedGene ~ overallMapRate, data=pd)

plot(totalAssignedGene ~ ERCCsumLogErr, data=pd)
plot(overallMapRate ~ ERCCsumLogErr, data=pd)
plot(mitoRate ~ ERCCsumLogErr, data=pd)

