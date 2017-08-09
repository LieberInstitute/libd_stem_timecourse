# qsub -tc 245 -l mf=45G,h_vmem=45G,h_stack=256M  -m e -M amanda.joy.price@gmail.com -t 1-57 IRFinder_add_57_stem_timecourse.sh

module load bedtools
module load samtools

export PATH=$PATH:/users/aprice26/biotools/STAR-2.5.2a/bin/Linux_x86_64

FILELIST=/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/IRFinder/shortIDs_57.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)

fastq_r1=/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/IRFinder/fastq_R1_57.txt
read1=$(awk "NR==$SGE_TASK_ID" $fastq_r1)

fastq_r2=/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/IRFinder/fastq_R2_57.txt
read2=$(awk "NR==$SGE_TASK_ID" $fastq_r2)

/users/aprice26/biotools/IRFinder-1.1.1/bin/IRFinder -d /dcl01/lieber/ajaffe/Amanda/libd_stem_timecourse/IRFinder/${ID} -r /dcl01/lieber/ajaffe/Amanda/NucVsCyt/Human-hg19-release75 $read1 $read2
