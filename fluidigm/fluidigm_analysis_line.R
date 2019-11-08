####
library(readxl)
library(lme4)
library(lmerTest)
library(RColorBrewer)
library(jaffelab)

dir.create("plots")

## phenotype
pd = read_excel("IPS_CODE_01252016.xlsx", sheet=1,skip=2)
pd = as.data.frame(pd[,1:6]) # keep first 6 columns

names(pd)[1] = "RunID"
pd = pd[!is.na(pd$RunID),]
pd$DONOR = ss(pd$RunID, "-|_", 1)
rownames(pd) = pd$RunID
pd$DONOR[pd$TYPE=="FIB"] = ss(pd$RunID[pd$TYPE=="FIB"], "-", 2)
pd$DONOR[pd$TYPE=="ES"] = "ES"

# make factors
pd$DONOR = factor(pd$DONOR)
pd$RPS = factor(pd$RPS, levels = c("LO", "HI"))
pd$DX = factor(pd$DX, levels = c("CNT", "SZ"))

# make anoter donor for modeling
pd$DONOR2 = as.character(pd$DONOR)
pd$DONOR2[pd$DONOR == "ES"] = pd$RunID[pd$DONOR == "ES"] 

# add colors
pd$pointCol = c(brewer.pal(9, "Set1"), brewer.pal(8,"Dark2"))[as.numeric(pd$DONOR)]
pd$pointShape = (23:21)[as.numeric(factor(pd$TYPE))]

# number of unique donors
length(unique(pd$DONOR[pd$TYPE %in% c("iPS","FIB")]))
length(unique(pd$DONOR[pd$TYPE == "iPS"]))
quantile(table(as.character(pd$DONOR[pd$TYPE == "iPS"])))

## expression
dat = as.data.frame(read_excel("FLUIDIGM_IPS_01252016.xlsx",sheet=1))

names(dat)[1] = "RunID"
rownames(dat) = dat$RunID
dat = dat[pd$RunID,]
dat$RunID = NULL
colnames(dat) = gsub(" ", "_", colnames(dat))
colnames(dat) = ss(colnames(dat), "\\.")

## normalize by HKs
hk = rowMeans(dat[,grep("GUSB|GAPDH", names(dat))])

# make into long file format
datLong = data.frame(ID = rep(rownames(dat), ncol(dat)),
	Gene = rep(colnames(dat), each=nrow(dat)), CT = unlist(dat))
datLong = datLong[order(datLong$Gene, datLong$ID),]
datLong$Rep = rep(1:2, times = nrow(datLong)/2)
datLong$CT[datLong$CT == 999] = 40 # make no expression = 40
datLong$HK = hk[match(datLong$ID, names(hk))]
datLong$deltaCT = 2^(datLong$HK - datLong$CT)

### mean mat for visualization
tmpList = split(datLong, factor(datLong$Gene,
	levels=unique(datLong$Gene)))
meanMat = sapply(tmpList, function(x) 
	tapply(x$deltaCT, x$ID, mean, na.rm=TRUE))
meanMat[is.nan(meanMat)] = NA
meanMat = meanMat[rownames(pd),]
meanMatRaw = sapply(tmpList, function(x) 
	tapply(x$CT, x$ID, mean, na.rm=TRUE))
meanMatRaw[is.nan(meanMatRaw)] = NA
meanMatRaw = meanMatRaw[rownames(pd),]

###############
# housekeeping check

## versus each other
par(mar=c(5,6,2,2), cex.axis=1.8, cex.lab=1.8)
plot(meanMatRaw[,"GUSB"], meanMatRaw[,"GAPDH"],
	pch = pd$pointShape, bg=pd$pointCol,cex=2,
	xlab="Mean GUSB", ylab="Mean GAPDH")
legend("topleft", levels(pd$DONOR),
	col=1:17,pch=15,nc=3, cex=1.6)
legend("bottomright", levels(factor(pd$TYPE)),	
	pch = 23:21,cex=2)

	
### very highly correlated
cor(meanMatRaw[,c("GUSB", "GAPDH")])^2
cor(meanMatRaw[-76,c("GUSB", "GAPDH")])^2

plot(hk, prcomp(dat[,grep("GUSB|GAPDH", names(dat))])$x[,1],
	pch = pd$pointShape, bg=pd$pointCol,cex=2,
	xlab="Mean HK", ylab="PC1: 99.2% Var Expl")
legend("topleft", levels(pd$DONOR),
	col=1:17,pch=15,nc=3, cex=1.6)
legend("bottomright", levels(factor(pd$TYPE)),	
	pch = 23:21,cex=2)

# by line
pdf("plots/rawHk_byDonor.pdf")
yl = range(meanMatRaw[,c("GUSB", "GAPDH")])
palette(c(brewer.pal(9, "Set1"), brewer.pal(8,"Dark2")))
par(mfrow = c(2,1), mar=c(3,6,4,2), cex.axis=1.8, cex.lab=1.8)
boxplot(meanMatRaw[,"GUSB"]~pd$DONOR,ylim = yl,
	las=3,xaxt="n", ylab="Mean GUSB",outline=FALSE)
points(meanMatRaw[,"GUSB"]~jitter(as.numeric(
	pd$DONOR),amount=0.15), pch=pd$pointShape, bg = pd$pointCol)
par(mar=c(8,6,0,2))
boxplot(meanMatRaw[,"GAPDH"]~pd$DONOR,
	ylim = yl,outline=FALSE,
	las=3,ylab="Mean GAPDH")
points(meanMatRaw[,"GAPDH"]~jitter(as.numeric(
	pd$DONOR),amount=0.15), pch=pd$pointShape, bg = pd$pointCol)
dev.off()

####################################
## drop the line w/ no HK expression
keepIndex=which(meanMatRaw[,"GAPDH"] < 30)
pd = pd[keepIndex,]
meanMat = meanMat[keepIndex,]
meanMatRaw = meanMatRaw[keepIndex,]
datLong = datLong[datLong$HK < 30,]

## versus each other
par(mar=c(5,6,2,2), cex.axis=1.8, cex.lab=1.8)
plot(meanMatRaw[,"GUSB"], meanMatRaw[,"GAPDH"],
	pch = pd$pointShape, bg=pd$pointCol,cex=2,
	xlab="Mean GUSB", ylab="Mean GAPDH")
legend("topleft", levels(pd$DONOR),
	col=1:17,pch=15,nc=3, cex=1.6)
legend("bottomright", levels(factor(pd$TYPE)),	
	pch = 23:21,cex=2)
	
pdf("plots/rawHk_byDonor_postdrop.pdf")
yl = range(meanMatRaw[,c("GUSB", "GAPDH")])
palette(c(brewer.pal(9, "Set1"), brewer.pal(8,"Dark2")))
par(mfrow = c(2,1), mar=c(3,6,4,2), cex.axis=1.8, cex.lab=1.8)
boxplot(meanMatRaw[,"GUSB"]~pd$DONOR,ylim = yl,
	las=3,xaxt="n", ylab="Mean GUSB",outline=FALSE)
points(meanMatRaw[,"GUSB"]~jitter(as.numeric(
	pd$DONOR),amount=0.15), pch=pd$pointShape, bg = pd$pointCol)
par(mar=c(8,6,0,2))
boxplot(meanMatRaw[,"GAPDH"]~pd$DONOR,
	ylim = yl,outline=FALSE,
	las=3,ylab="Mean GAPDH")
points(meanMatRaw[,"GAPDH"]~jitter(as.numeric(
	pd$DONOR),amount=0.15), pch=pd$pointShape, bg = pd$pointCol)
dev.off()

	
########################
### pca on mean mat ####
meanPca = prcomp(meanMat)
meanPcaRaw = prcomp(meanMatRaw)
meanPcs = as.data.frame(meanPca$x)
meanPcsRaw = as.data.frame(meanPcaRaw$x)
identical(rownames(pd), rownames(meanPcs)) # TRUE

pcaVars = getPcaVars(meanPca)
pcaVarsRaw = getPcaVars(meanPcaRaw)

pdf("plots/pca_byDonor.pdf", w=12)
par(mar=c(6,6,2,2))
sIndexes=split(1:nrow(pd), sort(pd$DONOR))
for(i in 1:ncol(meanPcs)) {
	plot(meanPcs[order(pd$DONOR),i], 
		pch =pd$pointShape[order(pd$DONOR)],
		bg = pd$pointCol[order(pd$DONOR)],
		cex=1.5, xaxt="n",xlab="DONOR",cex.axis=1.8,cex.lab=1.8,
		ylab=paste0("PC", i, ": ", pcaVars[i], "% Var Expl"))
	axis(1, at = sapply(sIndexes, mean), label = names(sIndexes),
		cex.axis=1.5, las=3)
	abline(v = sapply(sIndexes, max)+0.5, lty=2)
}
dev.off()

summary(lmer(meanPcs[,1] ~ TYPE + (1|DONOR2), data = pd)) # p=0.015

pdf("plots/pca_byDonor_raw.pdf", w=12)
palette(c(brewer.pal(9, "Set1"), brewer.pal(10,"Set3")))
par(mar=c(6,6,2,2))
for(i in 1:ncol(meanPcsRaw)) {
	plot(meanPcsRaw[order(pd$DONOR),i], 
		pch =pd$pointShape[order(pd$DONOR)],
		bg = pd$pointCol[order(pd$DONOR)],
		cex=1.5, xaxt="n",xlab="DONOR",cex.axis=1.8,cex.lab=1.8,
		ylab=paste0("PC", i, ": ", pcaVarsRaw[i], "% Var Expl"))
	axis(1, at = sapply(sIndexes, mean), label = names(sIndexes),
		cex.axis=1.5, las=3)
	abline(v = sapply(sIndexes, max)+0.5, lty=2)
}
dev.off()

summary(lmer(meanPcsRaw[,1] ~ TYPE + (1|DONOR2), 
	data = pd, subset=TYPE %in% c("iPS", "ES"))) # p=0.06

######################
# ANALYSIS FOR TYPE
#####################

## lmer
pd$TYPE = factor(pd$TYPE, levels = c("iPS", "FIB", "ES"))
coefListType = lapply(as.data.frame(meanMat),function(y) {
	try(summary(lmer(y ~ factor(TYPE) + 
		(1|DONOR),data = pd))$coef)
})

pvalMatType = t(sapply(coefListType, function(x) x[2:3,5]))
colnames(pvalMatType) = paste0(c("Fib", "ES"), "vsIPS_pval")
tMatType = t(sapply(coefListType, function(x) x[2:3,4]))
colnames(tMatType) = paste0(c("Fib", "ES"), "vsIPS_tstat")

sapply(c(0.05,0.01,0.001,1e-4), function(x) colSums(pvalMatType < x))
pvalMatType[pvalMatType[,1] > 0.05,]

outMatType = cbind(tMatType, pvalMatType)
outMatType = outMatType[,order(colnames(outMatType))]

### plots
pdf("plots/cellStateDiff.pdf",h=4.5,w=5)
par(mar=c(4,5,2,2))
for(i in 1:ncol(meanMat)) {
	boxplot(meanMat[,i] ~ pd$TYPE, las=3,
		xlab="",cex.axis=1.5,cex.lab=1.5,
		outline = FALSE, ylim = range(meanMat[,i]),
		cex.main=1.5, ylab="Norm Exprs (% of HK)",
		main = colnames(meanMat)[i])
	points(meanMat[,i] ~ jitter(
		as.numeric(pd$TYPE), amount=0.15), cex=1.5,
			pch = pd$pointShape, bg = pd$pointCol)
	text(2:3, y = max(meanMat[,i]), paste0("p=",
		signif(pvalMatType[i,],3)),cex=1.3)
}
dev.off()

#######################
### DONOR Analysis ####
#######################

ipsIndex=which(pd$TYPE=="iPS" & pd$DONOR != "62")
meanMatIps = meanMat[ipsIndex,]
pdIps = pd[ipsIndex,]
pdIps$DONOR = factor(as.character(pdIps$DONOR))

## difference in selection
pdIps$Donor_In_Analysis = pdIps$DONOR %in% c("165", "21", "3", "66", "90")

length(unique(pdIps$DONOR))
coefMatInc = t(sapply(as.data.frame(meanMatIps),function(y) {
	try(summary(lmer(y ~ factor(Donor_In_Analysis) +
		(1|DONOR),data = pdIps))$coef[2,])
}))
coefMatInc = as.data.frame(coefMatInc)
colnames(coefMatInc) = c("Diff", "SE", "df", "t", "p")
coefMatInc[coefMatInc$p < 0.01,]

## by subclone
pdIps$Clone_In_Analysis = pdIps$RunID %in% c( "165-B-2","165-B-3",
						"165-B-3X", "165-B-6X","21-B-3", "21-B-8", "21-B-9",
						"3-A-7", "3-A-8", "66-A-3", "66-A-9", "90-A-10", "90-A-5")
pdIpsDonor = pdIps[pdIps$Donor_In_Analysis,]
meanMapIpsDonor = meanMatIps[rownames(pdIpsDonor),]

coefMatClone = t(sapply(as.data.frame(meanMapIpsDonor),function(y) {
	try(summary(lmer(y ~ Clone_In_Analysis +
		(1|DONOR),data = pdIpsDonor))$coef[2,])
}))	
coefMatClone = as.data.frame(coefMatClone)
colnames(coefMatClone) = c("Diff", "SE", "df", "t", "p")
coefMatClone[coefMatClone$p < 0.01,]


## rerun PCA
meanPcaIps = prcomp(meanMatIps)
meanPcsIps = as.data.frame(meanPcaIps$x)
sIndexesIps=split(1:nrow(pdIps), sort(pdIps$DONOR))
pcaVarsIps = getPcaVars(meanPcaIps)

pdf("plots/pca_byDonor_iPSonly.pdf", w=10,h=6)
par(mar=c(6,6,2,2))
for(i in 1:ncol(meanPcsIps)) {
	plot(meanPcsIps[order(pdIps$DONOR),i], 
		pch =pdIps$pointShape[order(pdIps$DONOR)],
		bg = pdIps$pointCol[order(pdIps$DONOR)],
		cex=1.5, xaxt="n",xlab="DONOR",cex.axis=1.8,cex.lab=1.8,
		ylab=paste0("PC", i, ": ", pcaVarsIps[i], "% Var Expl"))
	axis(1, at = sapply(sIndexesIps, mean), label = names(sIndexesIps),
		cex.axis=1.5, las=3)
	abline(v = sapply(sIndexesIps, max)+0.5, lty=2)
}
dev.off()

### plots of PC1
ooPlot = order(abs(meanPcaIps$rot[,1]),decreasing=TRUE)[1:6]
pdf("plots/topPcRots_byDonor_iPSonly.pdf", w=7,h=6)
par(mar=c(6,6,3.5,2),cex.axis=1.8,cex.lab=1.8,cex.main=2.5)
for(i in ooPlot) {
	plot(meanMat[order(pdIps$DONOR),i], 
		pch =pdIps$pointShape[order(pdIps$DONOR)],
		bg = pdIps$pointCol[order(pdIps$DONOR)],
		cex=1.5, xaxt="n",xlab="DONOR",
		ylab="Norm Exprs (% of HK)", main = colnames(meanMat)[i])
	axis(1, at = sapply(sIndexesIps, mean), label = names(sIndexesIps),
		cex.axis=1.5, las=3)
	abline(v = sapply(sIndexesIps, max)+0.5, lty=2)
}
dev.off()

### hclust
dd = dist(meanPcsIps)
hc = hclust(dd)
myplclust(hc, lab.col=pdIps$pointCol,xlab="",sub="")

### variance components analysis ##
varCompAnalysis = apply(meanMatIps,2,function(y) {
	fit = lm(y ~ DONOR, data=pdIps)
	full = anova(fit)
	fullSS =full$"Sum Sq"
	signif(cbind(full,PctExp=fullSS/
		sum(fullSS)*100),3)
})
donorVarIps = t(sapply(varCompAnalysis, function(x) x[,6]))
colnames(donorVarIps) = gsub("factor(", "", 
	rownames(varCompAnalysis[[1]]), fixed=TRUE)
colnames(donorVarIps) = gsub(")", "",colnames(donorVarIps), fixed=TRUE)
colnames(donorVarIps) = paste0(colnames(donorVarIps), "_varExpl")
donorVarIps = as.data.frame(donorVarIps)
donorVarIps$donorPval = sapply(varCompAnalysis, function(x) x[1,5])

sapply(c(0.05,0.01,0.001,1e-4), function(x) sum(donorVarIps$donorPval < x))
cat(rownames(donorVarIps[donorVarIps$donorPval < 0.05,]))

rownames(donorVarIps[donorVarIps$donorPval < 0.05,]) %in% 
	rownames(pvalMat[pvalMat[,1] > 0.05,])

#### boxplots 
pdf("plots/donorLineDifference_iPSC_only.pdf", w=10,h=6)
par(mar=c(5,6,2,2))
for(i in 1:ncol(meanMatIps)) {
	boxplot(meanMatIps[,i] ~ pdIps$DONOR, las=3,
		xlab="",cex.axis=1.8,cex.lab=1.8,
		outline = FALSE, ylim = range(meanMatIps[,i]),
		cex.main=2, ylab="Norm Exprs (% of HK)",
		main = colnames(meanMatIps)[i])
	points(meanMatIps[,i] ~ jitter(
		as.numeric(pdIps$DONOR), amount=0.15),cex=1.7,
			pch = pdIps$pointShape, bg = pdIps$pointCol)
	legend("topright", paste0("p=", 
		signif(donorVarIps$donorPval[i], 3), "\n",
		donorVarIps$DONOR_varExpl[i],"% VE"),cex=1.8)
}
dev.off()

order(donorVarIps$donorPval)
donorVarIps[order(donorVarIps$donorPval)[1:6],]
write.csv(donorVarIps, file= "variance_explained_donor.csv")

#####################
### DX & RPS ########
####################

## add interaction
pdIps$dxRpsInt = as.numeric(pdIps$DX)*as.numeric(pdIps$RPS)
## fit model
coefListDxByRps = lapply(as.data.frame(meanMatIps),function(y) {
	summary(lmer(y ~ DX + RPS + # dxRpsInt + 
		(1|DONOR),data = pdIps))$coef
})

pvalMatDxByRps = t(sapply(coefListDxByRps, function(x) x[2:3,5]))
colnames(pvalMatDxByRps) = paste0(c("DX", "RPS"), "_pval")
tMatDxByRps = t(sapply(coefListDxByRps, function(x) x[2:3,4]))
colnames(tMatDxByRps) = paste0(c("DX", "RPS"), "_tstat")

outMatDxByRps = cbind(tMatDxByRps, pvalMatDxByRps)
outMatDxByRps = outMatDxByRps[,order(colnames(outMatDxByRps))]
outMatDxByRps = as.data.frame(outMatDxByRps)

#### boxplots 
pdIps$DxByRps = paste0(pdIps$DX, ":", pdIps$RPS)
pdIps$DxByRps = factor(pdIps$DxByRps,
	levels =c("CNT:LO", "CNT:HI",	"SZ:LO", "SZ:HI")) 

pdf("plots/dxByRpsDifference_iPSC_only.pdf", w=7,h=5)
par(mar=c(5,6,2,2))
for(i in 1:ncol(meanMatIps)) {
	boxplot(meanMatIps[,i] ~ pdIps$DxByRps, 
		xlab="",cex.axis=1.7,cex.lab=1.7,
		outline = FALSE, ylim = range(meanMatIps[,i]),
		cex.main=1.7, ylab="Norm Exprs (% of HK)",
		main = colnames(meanMatIps)[i])
	points(meanMatIps[,i] ~ jitter(
		as.numeric(pdIps$DxByRps), amount=0.15),cex=1.7,
			pch = pdIps$pointShape, bg = pdIps$pointCol)
	legend("topright", paste0(c("Dx", "RPS"), ": p=", 
		signif(pvalMatDxByRps[i,], 3)),cex=1.4)
}
dev.off()

## mean by donor
pdMeanIps = pdIps[!duplicated(pdIps$DONOR),]
dIndexesIps = split(1:nrow(pdIps), pdIps$DONOR)
meanByDonorIps = t(sapply(dIndexesIps, 
	function(ii) colMeans(meanMatIps[ii,])))
meanByDonorIps = meanByDonorIps[as.character(pdMeanIps$DONOR),]

pdf("plots/dxByRpsDifference_iPSC_only_means.pdf", w=7,h=5)
par(mar=c(5,6,2,2))
for(i in 1:ncol(meanByDonorIps)) {
	boxplot(meanByDonorIps[,i] ~ pdMeanIps$DxByRps, 
		xlab="",cex.axis=1.7,cex.lab=1.7,
		outline = FALSE, ylim = range(meanByDonorIps[,i]),
		cex.main=1.7, ylab="Norm Exprs (% of HK)",
		main = colnames(meanByDonorIps)[i])
	points(meanByDonorIps[,i] ~ jitter(
		as.numeric(pdMeanIps$DxByRps), amount=0.15),cex=1.7,
			pch = pdMeanIps$pointShape, bg = pdMeanIps$pointCol)
	legend("topright", paste0(c("Dx", "RPS"), ": p=", 
		signif(pvalMatDxByRps[i,], 3)),cex=1.4)
}
dev.off()

order(outMatDxByRps$DX_pval)[1:4]
order(outMatDxByRps$RPS_pval)[1:4]
