# !/bin/sh

# qsub -cwd -l mf=13G,h_vmem=20G,h_fsize=100G -t 2-20 download_reads.sh
# qsub -cwd -pe local 3 -l mf=5G,h_vmem=8G -t 2-41 download_reads.sh


# awk -F"\t" '{print $7}' ziller_SraRunTable.txt > srrs.txt
FILELIST=/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/ziller/srrs.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)
fastq-dump -O /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/ziller/FASTQ/ --gzip --split-files $ID 

