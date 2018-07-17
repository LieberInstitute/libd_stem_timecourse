#!/bin/bash
#$ -cwd
#$ -l mem_free=100G,h_vmem=100G,h_fsize=200G
#$ -N astro
#$ -o ./logs/astro_enrich_log.txt
#$ -e ./logs/astro_enrich_log.txt
#$ -m a
echo "**** Job starts ****"
date

Rscript /dcl01/lieber/ajaffe/lab/libd_stem_timecourse/tc_analysis/astrocytes/run_voom_neuron_astrocyte.R

echo "**** Job ends ****"
date
