##

pd = read.delim("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/ziller/ziller_SraRunTable.txt",
	as.is=TRUE)

readPath = "/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/ziller/FASTQ/"
man = data.frame(leftRead = paste0(readPath, pd$Run_s, "_1.fastq.gz"),
	leftMd5=0, rightRead = paste0(readPath, pd$Run_s, "_2.fastq.gz"),
	rightMd5= 0, SampleID = pd$Run_s,stringsAsFactors=FALSE)
write.table(man, file="samples.manifest",
	row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
	
###