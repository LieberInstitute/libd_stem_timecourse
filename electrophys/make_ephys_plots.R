##

ss = function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), "[", slot)

## read in data
library(readxl)
library(stringr)
library(rafalib)
library(RColorBrewer)
library(lmerTest)

## read in data
dat = read_excel("Neuron_analysis_Jaffe_clean.xlsx")
dat[dat %in% c("na","NA","")] = NA
dat = as.data.frame(dat)

dat$CellNumber = str_trim(dat$CellNumber)

## make numeric
for(i in 6:17) dat[,i] = as.numeric(dat[,i])

## add label
dat$Lab = paste0("Wk",dat$Week, "_", dat$Line)
dat$Donor = factor(ss(dat$Line, "A|B",1))

pal = c(colorRampPalette(brewer.pal(10,"Paired")[1:2])(5), 
		colorRampPalette(brewer.pal(10,"Paired")[3:4])(3), 
		brewer.pal(10,"Paired")[9:10], brewer.pal(10,"Paired")[5:8])
# 165B3 165B6  21B3  21B8  66A3  66A9 90A10  90A5
pal = pal[c(2,4, 6,7, 11:14)]
dat$lineCol = pal[as.numeric(factor(dat$Line))]

dat$donorWk= paste0(dat$Donor, ":", dat$Week,"wk")
colnames(dat)[6:7] = c("Capacitance", "MembraneResist")

## plots
pdf("../boxplots_of_activity_eeb.pdf",h=5,w=8,useDingbats=FALSE)
par(mar=c(7,6,3,2), cex.axis= 1.8, cex.lab=1.8, cex.main=1.6)
boxplot(dat$Capacitance ~ dat$donorWk,las=3,
	ylab = "Capacitance (pF)", outline=FALSE,
	ylim=range(dat$Capacitance, na.rm=TRUE))
points(dat$Capacitance ~ jitter(as.numeric(factor(dat$donorWk)),
	amount=0.15),pch = 21, bg = dat$lineCol)
boxplot(dat$MembraneResist ~ dat$donorWk,las=3,
	ylab = "Membr. Res. (MOhms)", outline=FALSE,
	ylim=c(0,2300))
points(dat$MembraneResist ~ jitter(as.numeric(factor(dat$donorWk)),
	amount=0.15),pch = 21, bg = dat$lineCol)
dev.off()

## capacitance, membrane resistance, peak Na current (ramp) pA
summCap = summary(lmer(log(Capacitance) ~ Week + (1|Donor) + (1|Line), data=dat))
summCap
summMem = summary(lmer(log(MembraneResist) ~ Week + (1|Donor) + (1|Line), data=dat))
summMem
#################
## transforms??##
mypar(3,4)
for(i in 6:17) hist(dat[,i],main = colnames(dat)[i])
for(i in 6:17) boxplot(dat[,i] ~ dat$Lab,
	main = colnames(dat)[i], las = 3)

datXform = as.matrix(dat[,6:17])
for(i in 1:ncol(datXform)) {
	x = datXform[,i]
	theSign = unique(sign(x))
	theSign = theSign[!is.na(theSign)]
	datXform[,i] = theSign*log(abs(x)+1)
}
datXform[!is.finite(datXform)] = NA

for(i in 1:12) hist(datXform[,i],main = colnames(dat)[i])
## looks better

## split into data and phenotype
pd = dat[,c(1:5,18,19)]
y  = as.data.frame(datXform)