#!/bin/bash
# Submit all chrom jobs for generation of gene-based windows.

windowSize=$1

########ADD IN CHROM BEDGRAPH GENRATION###########

#First, run window generation for all chromosomes in parallel.
for chrom in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
do
job_id_list+=$(sbatch --parsable --job-name=${chrom}_window gene_window_chr_submit.sh ${chrom} ${windowSize}):
done

echo "Submitted jobs:  ${job_id_list}"

#Second, run a script concatenating all chromosomes and applying pval correction
job_id_list=${job_id_list::-1} #removes superfluous semicolon at end of string
sbatch --dependency=afterok:${job_id_list} join_windows.sh ${windowSize}

echo "Submitted batch job ${join_job_id}"
