# qsub -tc 245 -l mf=50G,h_vmem=50G,h_stack=256M  -m e -M amanda.joy.price@gmail.com -t 1-106 IRFinder_libd_stem_timecourse_hg38_106.sh

module load bedtools
module load samtools

export PATH=$PATH:/users/aprice26/biotools/STAR-2.5.2a/bin/Linux_x86_64

FILELIST=/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/IRFinder/shortID_hg38_106.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)

fastq_r1=/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/IRFinder/fastq_R1_hg38_106.txt
read1=$(awk "NR==$SGE_TASK_ID" $fastq_r1)

fastq_r2=/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/IRFinder/fastq_R2_hg38_106.txt
read2=$(awk "NR==$SGE_TASK_ID" $fastq_r2)

echo ${ID} 

/users/aprice26/biotools/IRFinder-1.1.1/bin/IRFinder -d /dcl01/lieber/ajaffe/Amanda/libd_stem_timecourse/IRFinder/hg38/${ID} -r /users/aprice26/biotools/IRFinder-1.1.1/REF/Human-hg38-release81 $read1 $read2
