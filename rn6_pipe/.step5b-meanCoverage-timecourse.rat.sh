#!/bin/bash
#$ -cwd
#$ -l mem_free=100G,h_vmem=120G,h_fsize=100G
#$ -N step5b-meanCoverage-timecourse.rat
#$ -o ./logs/meanCoverage-timecourse.txt
#$ -e ./logs/meanCoverage-timecourse.txt
#$ -hold_jid pipeline_setup,step5-coverage-timecourse.rat
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

if [ ! -f "inferred_strandness_pattern.txt" ]
then
    echo "Missing the file inferred_strandness_pattern.txt"
    exit 1
fi

## Load required software
module load wiggletools/default
module load ucsctools

## Read the strandness information
STRANDRULE=$(cat inferred_strandness_pattern.txt)

if [ ${STRANDRULE} == "none" ]
then
    echo "*****************************"
    date
    echo "Processing unstranded data"
    echo "*****************************"
    ## Locate normalized BigWig files and concatenate them in a space separated list
    BIGWIGS=$(cat /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/samples.manifest | awk '{print "/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/Coverage/"$NF".bw"}' | paste -sd " ")
    
    ## Create mean of normalized bigwigs
    wiggletools write /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/Coverage/mean.wig mean ${BIGWIGS}
    wigToBigWig /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/Coverage/mean.wig /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/rn6.chrom.sizes.ensembl /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/Coverage/mean.bw
else
    for strand in Forward Reverse
    do
        echo "*****************************"
        date
        echo "Processing strand ${strand}"
        echo "*****************************"
        ## Locate normalized BigWig files and concatenate them in a space separated list
        BIGWIGS=$(cat /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/samples.manifest | awk -v strand="${strand}" '{print "/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/Coverage/"$NF"."strand".bw"}' | paste -sd " ")
    
        ## Create mean of normalized bigwigs
        wiggletools write /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/Coverage/mean.${strand}.wig mean ${BIGWIGS}
        wigToBigWig /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/Coverage/mean.${strand}.wig /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/rn6.chrom.sizes.ensembl /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/Coverage/mean.${strand}.bw
    done
fi

## Remove temp files
rm /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/rn6_pipe/Coverage/mean*.wig

echo "**** Job ends ****"
date
