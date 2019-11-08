#!/bin/bash
#$ -cwd
#$ -l mf=10G,h_vmem=10G,h_fsize=100G,h_stack=256M
#$ -N bulk_download_reads
#$ -o bulk/logs/download.$TASK_ID.txt
#$ -e bulk/logs/download.$TASK_ID.txt
#$ -t 1-456
#$ -tc 456

FILELIST=/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/volpato/volpato_SraRunTable_bulk.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST | cut -f9)

/users/ajaffe/software/sratoolkit.2.8.1-3-centos_linux64/bin/fastq-dump -O /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/volpato/bulk/FASTQ/ --gzip --split-files $ID
