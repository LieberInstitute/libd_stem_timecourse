###
library(SummarizedExperiment)
library(recount)
library(minfi)
library(GEOquery)
library(edgeR)
library(limma)
library(jaffelab)

# load data
load("preprocessed_human/rse_gene_tekin_human_n157.Rdata")
rse_gene$Species[rse_gene$Species == "Human and mouse"] = "Mix"

## get expression
yExprs = log2(getRPKM(rse_gene,length_var="Length")+1)

## figure out phenotypes
table(rse_gene$age_days, rse_gene$Species)

## figure out better fields
rse_gene$isSorted = grepl("sorted", rse_gene$strain_cell_type)

## geo phenotypes
geoData = getGEO("GSE111831")
geo_pheno = do.call("rbind", lapply(geoData, function(x) pData(x)[,c("geo_accession","title")]))

rse_gene$geo_title = as.character(geo_pheno$title[match(rse_gene$Sample_Name, 
					geo_pheno$geo_accession)])
rse_gene$geo_label = ss(rse_gene$geo_title, "_rep")

###############################
## check composition first ####
###############################

# project
load("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/deconvolution/cell_type/singleCell_iPSC_quake_coefEsts_calibration_Zscale.rda")
yExprs_Z = scale(yExprs[rownames(coefEsts),])
propEsts = as.data.frame(minfi:::projectCellType(yExprs_Z,coefEsts))

boxplot(propEsts$Neuron ~ rse_gene$age_days)
boxplot(propEsts$Neuron ~ rse_gene$culture_condition)
plot(propEsts$Neuron ~ rse_gene$age_days)

par(mar=c(5,30,1,1))
boxplot(propEsts$iPSC ~ rse_gene$geo_label,horizontal=TRUE,las=1)
boxplot(propEsts$NPC ~ rse_gene$geo_label,horizontal=TRUE,las=1)
boxplot(propEsts$Fetal_quiescent ~ rse_gene$geo_label,horizontal=TRUE,las=1)
boxplot(propEsts$Neuron ~ rse_gene$geo_label,horizontal=TRUE,las=1)

boxplot(propEsts$Astrocytes ~ rse_gene$geo_label,horizontal=TRUE,las=1)
boxplot(propEsts$Oligodendrocytes ~ rse_gene$geo_label,horizontal=TRUE,las=1)
boxplot(propEsts$Endothelial ~ rse_gene$geo_label,horizontal=TRUE,las=1)

######################
## plots for paper  ##
######################

groups = c("2D_1wk_iN","2D_1wk_iN_ma", "2D_5wk_iN", "2D_5wk_iN_ma",
	"iN_5wk_3D_100ul_7.36mg/ml_Mtrgl_10mpml",
	"iN_from_5wk_3D_iN_astrocytic-cell_100ul_7.36mg/ml_Mtrgl_20mpml",
	"iN_from_5wk_3D_iN_human-primary-astro_100ul_7.36mg/ml_Mtrgl_20mpml",
	
################################
## astrocytes ##################
################################

astro = read.csv("tc_analysis/astrocytes/astrocytes_day77_top_effects.csv")