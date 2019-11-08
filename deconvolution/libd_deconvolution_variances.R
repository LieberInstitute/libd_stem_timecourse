###
library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(genefilter)
library(RColorBrewer)
library(GEOquery)

MAINDIR = "/dcl01/lieber/ajaffe/lab/libd_stem_timecourse"



###############################################################
############## brain stage ####################################

pal = brewer.pal(8,"Dark2")

### LIBD stemcell
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n146.rda"))
## don't use technical replicates
bioInd = which(colData(rse_gene)$Class == "Naked genomes")
rse_gene = rse_gene[,bioInd]  # n=128
colData(rse_gene) = colData(rse_gene)
colData(rse_gene)$DAY = paste0("day_",colData(rse_gene)$DAY)
colData(rse_gene)$experiment = "lieber"
rse_geneTC = rse_gene

## clean CONDITION
rse_geneTC$COND = rse_geneTC$CONDITION
rse_geneTC$COND[grep("ACC_DORSAL",rse_geneTC$COND)] = "ACC_DORSAL"
rse_geneTC$COND[grep("PLUS",rse_geneTC$COND)] = "NEURONS+ASTROS"
rse_geneTC$COND = as.factor(rse_geneTC$COND)
rse_geneTC$COND = factor(rse_geneTC$COND,levels(rse_geneTC$COND)[c(5,1,4,6,2,3)])

geneRpkmTC = recount::getRPKM(rse_geneTC, "Length")

load("brain_stage/iPSC_BrainSpan_coefEsts_calibration_Zscore.rda")

yExprsTC = log2(geneRpkmTC+1)
yExprsTC_Z = scale(yExprsTC[rownames(coefEsts),])

### fit model
stemPropEsts = minfi:::projectCellType(yExprsTC_Z,coefEsts)
stemPropEsts_Prop = prop.table(stemPropEsts,1)

### plots
gIndexes_tc = splitit(rse_geneTC$COND)
stemPropEsts_groupMeans_tc = sapply(gIndexes_tc, 
	function(ii) colMeans(stemPropEsts[ii,]))
stemPropEsts_groupSEs_tc = sapply(gIndexes_tc, 
	function(ii) apply(stemPropEsts[ii,],2,function(x) sd(x)/sqrt(length(x))))

# for the paper
round(100*stemPropEsts_groupMeans_tc,2)
round(100*stemPropEsts_groupSEs_tc,2)


### Tekin

# load data
load("../otherData/Tekin_nature/preprocessed_human/rse_gene_tekin_human_n157.Rdata")
rse_gene_tekin = rse_gene

## geo phenotypes
geoData = getGEO("GSE111831")
geo_pheno = do.call("rbind", lapply(geoData, function(x) pData(x)[,c("geo_accession","title")]))
rse_gene_tekin$geo_title = as.character(geo_pheno$title[match(rse_gene_tekin$Sample_Name, 
					geo_pheno$geo_accession)])
rse_gene_tekin$geo_label = ss(rse_gene_tekin$geo_title, "_rep")


type = c("human neuronal cells sorted from their 3D co-culture with human primary astrocytes",
		"human neuronal cells sorted from their 3D co-culture with differentiated human astrocytic cells",
		"human neuronal cells sorted from 3D only human neuronal cell culture",
		"human neuronal cells and mouse astrocytes",
		"human neuronal cells")
keepInd = which(rse_gene_tekin$strain_cell_type %in% type &
				rse_gene_tekin$culture_condition %in% c("2D","3D") &
				rse_gene_tekin$age %in% c("1wk","5wk"))			
rse_gene_tekin = rse_gene_tekin[,keepInd]

rse_gene_tekin$cells = ifelse(grepl("astro", rse_gene_tekin$strain_cell_type), "neur+ast", "neur")
rse_gene_tekin$SampleLabel =paste0(rse_gene_tekin$age, "_", rse_gene_tekin$culture_condition, "_", rse_gene_tekin$cells)

dropInd = which(rse_gene_tekin$SampleLabel == "5wk_3D_neur+ast" & 
			rse_gene_tekin$strain_cell_type == "human neuronal cells and mouse astrocytes")
rse_gene_tekin = rse_gene_tekin[,-dropInd]

yExprs_tekin = log2(recount::getRPKM(rse_gene_tekin,length_var="Length")+1)


load("brain_stage/iPSC_BrainSpan_coefEsts_calibration_Zscore.rda")

yExprsTEK_Z = scale(yExprs_tekin[rownames(coefEsts),])

### fit model
stemPropEsts = minfi:::projectCellType(yExprsTEK_Z,coefEsts)
stemPropEsts_Prop = prop.table(stemPropEsts,1)

### plots
gIndexes_TEK = splitit(rse_gene_tekin$SampleLabel)
stemPropEsts_groupMeans_TEK = sapply(gIndexes_TEK, 
	function(ii) colMeans(stemPropEsts[ii,]))
stemPropEsts_groupSEs_TEK = sapply(gIndexes_TEK, 
	function(ii) apply(stemPropEsts[ii,],2,function(x) sd(x)/sqrt(length(x))))

round(100*stemPropEsts_groupMeans_TEK,2)
round(100*stemPropEsts_groupSEs_TEK,2)


### Bardy

pd = read.csv("../otherData/Bardy/metadata_56cells.csv",as.is=TRUE)

load("../otherData/Bardy/rse_gene_Bardy_singleCell_n56.Rdata")
rse_gene_bardy = rse_gene
rse_gene_bardy = rse_gene_bardy[,pd$cell.number]
colData(rse_gene_bardy) = cbind(colData(rse_gene_bardy), pd)
yExprs_bardy = log2(recount::getRPKM(rse_gene_bardy,length_var="Length")+1)


load("brain_stage/iPSC_BrainSpan_coefEsts_calibration_Zscore.rda")

yExprsBardy_Z = scale(yExprs_bardy[rownames(coefEsts),])

### fit model
stemPropEsts = minfi:::projectCellType(yExprsBardy_Z,coefEsts)
stemPropEsts_Prop = prop.table(stemPropEsts,1)

### plots
gIndexes_Bardy = splitit(rse_gene_bardy$FTC_simple)
stemPropEsts_groupMeans_Bardy = sapply(gIndexes_Bardy, 
	function(ii) colMeans(stemPropEsts[ii,]))
stemPropEsts_groupSEs_Bardy = sapply(gIndexes_Bardy, 
	function(ii) apply(stemPropEsts[ii,],2,function(x) sd(x)/sqrt(length(x))))

round(100*stemPropEsts_groupMeans_Bardy,2)
round(100*stemPropEsts_groupSEs_Bardy,2)






###############################################################
############## cell type ######################################

pal = brewer.pal(10, "Set3")
pal[2] = "#8B795E" # darken NPCs


### LIBD stemcell
load(file.path(MAINDIR,"data/libd_stemcell_timecourse_rseGene_n146.rda"))
## don't use technical replicates
bioInd = which(colData(rse_gene)$Class == "Naked genomes")
rse_gene = rse_gene[,bioInd]  # n=128
colData(rse_gene) = colData(rse_gene)
colData(rse_gene)$DAY = paste0("day_",colData(rse_gene)$DAY)
colData(rse_gene)$experiment = "lieber"
rse_geneTC = rse_gene

## clean CONDITION
rse_geneTC$COND = rse_geneTC$CONDITION
rse_geneTC$COND[grep("ACC_DORSAL",rse_geneTC$COND)] = "ACC_DORSAL"
rse_geneTC$COND[grep("PLUS",rse_geneTC$COND)] = "NEURONS+ASTROS"
rse_geneTC$COND = as.factor(rse_geneTC$COND)
rse_geneTC$COND = factor(rse_geneTC$COND,levels(rse_geneTC$COND)[c(5,1,4,6,2,3)])

geneRpkmTC = recount::getRPKM(rse_geneTC, "Length")

load("cell_type/singleCell_iPSC_quake_coefEsts_calibration_Zscale.rda")

yExprsTC = log2(geneRpkmTC+1)
yExprsTC_Z = scale(yExprsTC[rownames(coefEsts),])

### fit model
stemPropEsts = minfi:::projectCellType(yExprsTC_Z,coefEsts)
stemPropEsts_Prop = prop.table(stemPropEsts,1)

### plots
gIndexes_tc = splitit(rse_geneTC$COND)
stemPropEsts_groupMeans_tc_celltype = sapply(gIndexes_tc, 
	function(ii) colMeans(stemPropEsts[ii,]))
stemPropEsts_groupSEs_tc_celltype = sapply(gIndexes_tc, 
	function(ii) apply(stemPropEsts[ii,],2,function(x) sd(x)/sqrt(length(x))))

# for the paper
round(100*stemPropEsts_groupMeans_tc_celltype,2)
round(100*stemPropEsts_groupSEs_tc_celltype,2)


### Tekin

# load data
load("../otherData/Tekin_nature/preprocessed_human/rse_gene_tekin_human_n157.Rdata")
rse_gene_tekin = rse_gene

## geo phenotypes
geoData = getGEO("GSE111831")
geo_pheno = do.call("rbind", lapply(geoData, function(x) pData(x)[,c("geo_accession","title")]))
rse_gene_tekin$geo_title = as.character(geo_pheno$title[match(rse_gene_tekin$Sample_Name, 
					geo_pheno$geo_accession)])
rse_gene_tekin$geo_label = ss(rse_gene_tekin$geo_title, "_rep")


type = c("human neuronal cells sorted from their 3D co-culture with human primary astrocytes",
		"human neuronal cells sorted from their 3D co-culture with differentiated human astrocytic cells",
		"human neuronal cells sorted from 3D only human neuronal cell culture",
		"human neuronal cells and mouse astrocytes",
		"human neuronal cells")
keepInd = which(rse_gene_tekin$strain_cell_type %in% type &
				rse_gene_tekin$culture_condition %in% c("2D","3D") &
				rse_gene_tekin$age %in% c("1wk","5wk"))			
rse_gene_tekin = rse_gene_tekin[,keepInd]

rse_gene_tekin$cells = ifelse(grepl("astro", rse_gene_tekin$strain_cell_type), "neur+ast", "neur")
rse_gene_tekin$SampleLabel =paste0(rse_gene_tekin$age, "_", rse_gene_tekin$culture_condition, "_", rse_gene_tekin$cells)

dropInd = which(rse_gene_tekin$SampleLabel == "5wk_3D_neur+ast" & 
			rse_gene_tekin$strain_cell_type == "human neuronal cells and mouse astrocytes")
rse_gene_tekin = rse_gene_tekin[,-dropInd]

yExprs_tekin = log2(recount::getRPKM(rse_gene_tekin,length_var="Length")+1)


load("cell_type/singleCell_iPSC_quake_coefEsts_calibration_Zscale.rda")

yExprsTEK_Z = scale(yExprs_tekin[rownames(coefEsts),])

### fit model
stemPropEsts = minfi:::projectCellType(yExprsTEK_Z,coefEsts)
stemPropEsts_Prop = prop.table(stemPropEsts,1)

### plots
gIndexes_TEK = splitit(rse_gene_tekin$SampleLabel)
stemPropEsts_groupMeans_TEK_celltype = sapply(gIndexes_TEK, 
	function(ii) colMeans(stemPropEsts[ii,]))
stemPropEsts_groupSEs_TEK_celltype = sapply(gIndexes_TEK, 
	function(ii) apply(stemPropEsts[ii,],2,function(x) sd(x)/sqrt(length(x))))

round(100*stemPropEsts_groupMeans_TEK_celltype,2)
round(100*stemPropEsts_groupSEs_TEK_celltype,2)


### Bardy

pd = read.csv("../otherData/Bardy/metadata_56cells.csv",as.is=TRUE)

load("../otherData/Bardy/rse_gene_Bardy_singleCell_n56.Rdata")
rse_gene_bardy = rse_gene
rse_gene_bardy = rse_gene_bardy[,pd$cell.number]
colData(rse_gene_bardy) = cbind(colData(rse_gene_bardy), pd)
yExprs_bardy = log2(recount::getRPKM(rse_gene_bardy,length_var="Length")+1)


load("cell_type/singleCell_iPSC_quake_coefEsts_calibration_Zscale.rda")

yExprsBardy_Z = scale(yExprs_bardy[rownames(coefEsts),])

### fit model
stemPropEsts = minfi:::projectCellType(yExprsBardy_Z,coefEsts)
stemPropEsts_Prop = prop.table(stemPropEsts,1)

### plots
gIndexes_Bardy = splitit(rse_gene_bardy$FTC_simple)
stemPropEsts_groupMeans_Bardy_celltype = sapply(gIndexes_Bardy, 
	function(ii) colMeans(stemPropEsts[ii,]))
stemPropEsts_groupSEs_Bardy_celltype = sapply(gIndexes_Bardy, 
	function(ii) apply(stemPropEsts[ii,],2,function(x) sd(x)/sqrt(length(x))))

round(100*stemPropEsts_groupMeans_Bardy_celltype,2)
round(100*stemPropEsts_groupSEs_Bardy_celltype,2)












######## Results ##########


## Cell Type

> round(100*stemPropEsts_groupMeans_tc_celltype,2)
                  RENEW ACC_DORSAL   NPC ROSETTE NEURONS_ALONE NEURONS+ASTROS
iPSC              97.61      58.97 25.46   19.23          4.41           3.17
NPC                0.79      40.59 71.25   71.34         53.53          38.66
Fetal_replicating  0.00       0.15  1.63    3.70          0.00           0.00
Fetal_quiescent    2.22       1.46  1.82    7.87         21.31          32.98
OPC                0.00       0.00  0.00    0.00          0.00           0.00
Neurons            0.19       1.67  0.06    2.77         23.40          48.89
Astrocytes         0.01       0.16  0.38    0.14         15.96           1.25
Oligodendrocytes   0.04       0.00  0.00    0.00          0.00           0.00
Microglia          0.00       0.00  0.00    0.00          0.00           0.00
Endothelial        2.89       6.20 15.89   13.57         27.07          11.28
> round(100*stemPropEsts_groupSEs_tc_celltype,2)
                  RENEW ACC_DORSAL  NPC ROSETTE NEURONS_ALONE NEURONS+ASTROS
iPSC               0.50       3.76 0.84    0.96          1.58           0.70
NPC                0.40       3.84 1.55    1.32          2.56           0.63
Fetal_replicating  0.00       0.08 0.87    1.26          0.00           0.00
Fetal_quiescent    0.40       0.22 1.22    2.44          2.73           1.04
OPC                0.00       0.00 0.00    0.00          0.00           0.00
Neurons            0.10       0.23 0.04    1.10          1.07           1.41
Astrocytes         0.01       0.06 0.19    0.14          3.55           0.46
Oligodendrocytes   0.03       0.00 0.00    0.00          0.00           0.00
Microglia          0.00       0.00 0.00    0.00          0.00           0.00
Endothelial        0.46       0.35 1.71    3.25          6.25           2.32

> round(100*stemPropEsts_groupMeans_TEK_celltype,2)
                  1wk_2D_neur 1wk_2D_neur+ast 1wk_3D_neur 1wk_3D_neur+ast 5wk_2D_neur 5wk_2D_neur+ast 5wk_3D_neur 5wk_3D_neur+ast
iPSC                    17.19           14.41       18.01           15.02       10.77            7.86        9.68            7.50
NPC                     42.87           46.99       43.30           48.37       24.83           33.86       29.22           31.30
Fetal_replicating        0.00            0.00        0.00            0.00        0.00            0.00        0.00            0.00
Fetal_quiescent         34.01           27.41       28.95           23.86       41.01           32.56       40.62           32.43
OPC                      0.00            0.00        0.00            0.00        0.00            0.00        0.00            0.45
Neurons                 32.29           33.58       36.74           34.76       51.96           49.86       49.17           56.08
Astrocytes               0.00            0.00        0.00            0.00        0.00            0.00        0.00            0.00
Oligodendrocytes         0.00            0.00        0.00            0.00        0.00            0.00        0.30            1.18
Microglia                0.07            0.00        2.25            0.85        0.88            2.58        2.88            7.45
Endothelial              0.00            5.15        0.00            5.77        0.00            7.11        0.00            4.47
> round(100*stemPropEsts_groupSEs_TEK_celltype,2)
                  1wk_2D_neur 1wk_2D_neur+ast 1wk_3D_neur 1wk_3D_neur+ast 5wk_2D_neur 5wk_2D_neur+ast 5wk_3D_neur 5wk_3D_neur+ast
iPSC                     0.33            1.10        0.75            0.52        0.43            0.72        0.52            0.42
NPC                      0.31            1.34        0.73            0.51        0.28            0.19        0.90            1.36
Fetal_replicating        0.00            0.00        0.00            0.00        0.00            0.00        0.00            0.00
Fetal_quiescent          0.80            0.72        0.92            0.86        0.57            0.71        1.54            0.56
OPC                      0.00            0.00        0.00            0.00        0.00            0.00        0.00            0.28
Neurons                  0.71            0.60        0.97            0.84        0.65            0.33        1.03            1.13
Astrocytes               0.00            0.00        0.00            0.00        0.00            0.00        0.00            0.00
Oligodendrocytes         0.00            0.00        0.00            0.00        0.00            0.00        0.20            0.76
Microglia                0.04            0.00        0.15            0.25        0.35            0.55        0.60            0.80
Endothelial              0.00            0.75        0.00            1.06        0.00            0.57        0.00            1.75

> round(100*stemPropEsts_groupMeans_Bardy_celltype,2)
                      0     1     2     3     4     5
iPSC               1.74  4.38  6.28  3.87  6.13  6.39
NPC               59.40 41.87 35.77 33.18 30.53 38.90
Fetal_replicating  0.32  0.37  1.29  0.94  0.00  0.16
Fetal_quiescent    5.83 36.55 39.77 34.04 41.14 33.94
OPC                6.13  0.00  1.50  0.85  0.00  1.32
Neurons            3.35 15.07 11.71 23.72 30.78 30.99
Astrocytes        23.68  1.08  3.88  8.56  0.18  2.31
Oligodendrocytes   1.15  0.00  0.46  0.28  0.65  0.73
Microglia          2.44  3.39  5.88  2.82  4.91  5.88
Endothelial       20.01  7.22  9.62 10.96  3.86  5.56
> round(100*stemPropEsts_groupSEs_Bardy_celltype,2)
                     0    1    2    3    4    5
iPSC              1.29 1.43 1.11 1.54 1.72 1.04
NPC               2.18 8.35 4.81 6.70 3.17 2.80
Fetal_replicating 0.32 0.23 0.65 0.94 0.00 0.12
Fetal_quiescent   2.86 6.48 5.49 5.23 2.48 2.70
OPC               2.30 0.00 1.50 0.85 0.00 0.61
Neurons           1.82 5.05 3.88 6.18 2.58 2.31
Astrocytes        3.42 0.66 3.13 4.94 0.18 0.81
Oligodendrocytes  0.83 0.00 0.46 0.28 0.65 0.42
Microglia         1.41 2.06 2.00 1.69 1.58 0.71
Endothelial       2.01 2.73 2.09 3.77 1.94 1.20






## Brain Stage



> round(100*stemPropEsts_groupMeans_tc,2)
               RENEW ACC_DORSAL   NPC ROSETTE NEURONS_ALONE NEURONS+ASTROS
iPSC           96.14      75.06 56.27   49.15         23.75          27.52
NCX_EarlyFetal  0.46      20.58 33.27   33.27         16.73          29.35
NCX_MidFetal    0.00       0.12  0.00    9.63         33.47          29.88
NCX_LateFetal   0.00       1.18 10.50   10.91         10.10           1.99
NCX_Infant      0.00       1.14  0.25    0.39         16.37           3.71
NCX_Child       0.00       0.00  0.00    0.00          0.00           0.00
NCX_Teens       0.00       0.08  0.00    0.00          0.00           0.79
NCX_Adult       0.00       0.01  0.00    0.00          2.74          22.23
> round(100*stemPropEsts_groupSEs_tc,2)
               RENEW ACC_DORSAL  NPC ROSETTE NEURONS_ALONE NEURONS+ASTROS
iPSC            0.19       2.21 0.99    1.66          2.30           0.70
NCX_EarlyFetal  0.19       2.14 2.52    2.29          4.97           1.80
NCX_MidFetal    0.00       0.10 0.00    2.99          7.84           1.78
NCX_LateFetal   0.00       0.25 1.40    2.40          6.00           1.14
NCX_Infant      0.00       0.29 0.25    0.39          3.81           0.80
NCX_Child       0.00       0.00 0.00    0.00          0.00           0.00
NCX_Teens       0.00       0.04 0.00    0.00          0.00           0.79
NCX_Adult       0.00       0.01 0.00    0.00          1.37           1.06

> round(100*stemPropEsts_groupMeans_TEK,2)
               1wk_2D_neur 1wk_2D_neur+ast 1wk_3D_neur 1wk_3D_neur+ast 5wk_2D_neur 5wk_2D_neur+ast 5wk_3D_neur 5wk_3D_neur+ast
iPSC                 37.71           36.36       40.37           40.86       27.60           29.03       30.20           29.34
NCX_EarlyFetal       39.36           34.64       35.62           35.81       22.68           18.48       22.41           25.91
NCX_MidFetal         12.76           16.91       12.64           10.87       34.34           31.71       32.46           23.96
NCX_LateFetal         0.00            0.00        0.00            0.00        0.00            0.00        0.00            2.29
NCX_Infant           19.72           17.85       20.00           24.37       11.39           20.85       10.23            4.34
NCX_Child             0.00            0.00        0.00            0.00        0.00            0.00        0.00            0.00
NCX_Teens             4.19            8.15        2.35            8.35        0.00            1.99        0.20            0.00
NCX_Adult             3.43            4.59        6.36            0.00       21.25           18.21       22.61           28.34
> round(100*stemPropEsts_groupSEs_TEK,2)
               1wk_2D_neur 1wk_2D_neur+ast 1wk_3D_neur 1wk_3D_neur+ast 5wk_2D_neur 5wk_2D_neur+ast 5wk_3D_neur 5wk_3D_neur+ast
iPSC                  0.30            0.16        0.85            0.41        0.33            0.21        0.39            1.05
NCX_EarlyFetal        3.26            2.17        1.34            1.25        1.80            2.47        1.80            2.68
NCX_MidFetal          3.45            2.72        1.29            1.39        2.06            2.75        1.96            2.61
NCX_LateFetal         0.00            0.00        0.00            0.00        0.00            0.00        0.00            1.45
NCX_Infant            1.51            1.45        1.17            1.09        1.21            1.49        1.20            2.67
NCX_Child             0.00            0.00        0.00            0.00        0.00            0.00        0.00            0.00
NCX_Teens             3.04            4.08        1.60            1.01        0.00            1.99        0.20            0.00
NCX_Adult             1.91            4.59        1.72            0.00        0.46            2.34        1.31            2.24

> round(100*stemPropEsts_groupMeans_Bardy,2)
                   0     1     2     3     4     5
iPSC           30.61 33.32 34.20 32.40 35.03 34.58
NCX_EarlyFetal 12.55 35.22 30.41 25.11 36.16 34.72
NCX_MidFetal   17.85  4.80  9.55  9.34  6.45  6.96
NCX_LateFetal   3.63  0.18  0.00  1.80  0.00  0.28
NCX_Infant     22.41 20.81 21.98 23.12 20.89 25.39
NCX_Child      15.36 11.42 12.94 15.37  9.39  7.19
NCX_Teens       0.00  1.66  1.94  2.01  4.40  3.88
NCX_Adult       2.20  2.13  0.97  1.86  0.88  1.91
> round(100*stemPropEsts_groupSEs_Bardy,2)
                  0    1    2    3    4    5
iPSC           1.60 2.54 1.31 1.29 1.13 0.61
NCX_EarlyFetal 5.85 4.12 5.21 4.68 3.44 1.97
NCX_MidFetal   5.97 2.22 6.27 4.71 3.58 1.71
NCX_LateFetal  2.44 0.18 0.00 1.80 0.00 0.28
NCX_Infant     3.00 5.36 4.54 3.71 3.92 1.63
NCX_Child      4.67 3.24 3.64 4.44 4.69 1.41
NCX_Teens      0.00 1.66 1.42 2.01 1.50 1.18
NCX_Adult      2.20 2.13 0.97 1.86 0.60 0.88












