####
library(jaffelab)

## reads
readPath = "/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/otherData/Bardy/FASTQ"
fq = list.files(readPath, full.names=TRUE)
names(fq) = list.files(readPath)
names(fq) = ss(ss(names(fq),"_"), "-")

fqList = split(fq, names(fq))
fqMat = do.call("rbind", fqList)

manifest = data.frame(leftRead = fqMat[,1], leftMd5 = 0,
	rightRead = fqMat[,2], rightMd5 = 0, 
	SampleID = rownames(fqMat), stringsAsFactors=FALSE)
write.table(manifest, file = "samples.manifest",
	row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)