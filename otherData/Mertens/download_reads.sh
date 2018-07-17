#!/bin/bash
#$ -cwd
#$ -l mf=10G,h_vmem=10G,h_fsize=100G,h_stack=256M
#$ -N download_reads
#$ -o ./logs/download.$TASK_ID.txt
#$ -e ./logs/download.$TASK_ID.txt
#$ -t 2-19
#$ -tc 60

# https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP135684

FILELIST=/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/otherData/Mertens/Mertens_SraRunTable.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST | cut -f8)

/users/ajaffe/software/sratoolkit.2.8.1-3-centos_linux64/bin/fastq-dump -O /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/otherData/Mertens/FASTQ/ --gzip --split-files $ID
