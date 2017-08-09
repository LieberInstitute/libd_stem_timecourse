# qsub -tc 245 -l mf=45G,h_vmem=45G,h_stack=256M  -m e -M amanda.joy.price@gmail.com -t 1-112 IRFinder_libd_stem_timecourse.sh

module load bedtools
module load samtools

export PATH=$PATH:/users/aprice26/biotools/STAR-2.5.2a/bin/Linux_x86_64

FILELIST=/dcl01/ajaffe/data/Nina/AZ_Pilot/rids_dataNew_112.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)
read1=/dcl01/ajaffe/data/Nina/AZ_Pilot/data_new/${ID}_*_R1_001.fastq.gz
read2=/dcl01/ajaffe/data/Nina/AZ_Pilot/data_new/${ID}_*_R2_001.fastq.gz

/users/aprice26/biotools/IRFinder-1.1.1/bin/IRFinder -d /dcl01/lieber/ajaffe/Amanda/libd_stem_timecourse/IRFinder/${ID} -r /dcl01/lieber/ajaffe/Amanda/NucVsCyt/Human-hg19-release75 $read1 $read2