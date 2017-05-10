# qsub -tc 245 -l mf=60G,h_vmem=60G,h_stack=256M  -m e -M amanda.joy.price@gmail.com -t 1-11 IRFinder_repeat_11_more_memory.sh

module load bedtools
module load samtools

export PATH=$PATH:/users/aprice26/biotools/STAR-2.5.2a/bin/Linux_x86_64

FILELIST=/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/IRFinder/repeat_11_more_memory.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)
read1=/dcl01/ajaffe/data/Nina/AZ_Pilot/data_new/${ID}_*_R1_001.fastq.gz
read2=/dcl01/ajaffe/data/Nina/AZ_Pilot/data_new/${ID}_*_R2_001.fastq.gz

/users/aprice26/biotools/IRFinder-1.1.1/bin/IRFinder -d /dcl01/lieber/ajaffe/Amanda/libd_stem_timecourse/IRFinder/${ID} -r /dcl01/lieber/ajaffe/Amanda/NucVsCyt/Human-hg19-release75 $read1 $read2
