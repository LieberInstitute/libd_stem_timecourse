# qsub -tc 245 -l mf=60G,h_vmem=60G,h_stack=256M  -m e -M amanda.joy.price@gmail.com -t 1 repeat_R16-908_H7V3TBBXX.sh

module load bedtools
module load samtools

export PATH=$PATH:/users/aprice26/biotools/STAR-2.5.2a/bin/Linux_x86_64

rm -r /dcl01/lieber/ajaffe/Amanda/libd_stem_timecourse/IRFinder/R16-908_H7V3TBBXX
/users/aprice26/biotools/IRFinder-1.1.1/bin/IRFinder -d /dcl01/lieber/ajaffe/Amanda/libd_stem_timecourse/IRFinder/R16-908_H7V3TBBXX -r /dcl01/lieber/ajaffe/Amanda/NucVsCyt/Human-hg19-release75 /dcl01/ajaffe/data/Nina/AZ_Pilot/data/R16-908_H7V3TBBXX_S33_L007_R1_001.fastq.gz /dcl01/ajaffe/data/Nina/AZ_Pilot/data/R16-908_H7V3TBBXX_S33_L007_R2_001.fastq.gz
