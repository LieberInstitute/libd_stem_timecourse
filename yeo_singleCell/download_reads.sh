#!/bin/bash
#$ -cwd
#$ -l mf=10G,h_vmem=10G,h_fsize=100G,h_stack=256M
#$ -N download_reads
#$ -o ./logs/download.$TASK_ID.txt
#$ -e ./logs/download.$TASK_ID.txt
#$ -t 2-215
#$ -tc 60

FILELIST=/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/yeo_singleCell/yeo_SraRunTable.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST | cut -f8)

/users/ajaffe/software/sratoolkit.2.8.1-3-centos_linux64/bin/fastq-dump -O /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/yeo_singleCell/FASTQ/ --gzip --split-files $ID
