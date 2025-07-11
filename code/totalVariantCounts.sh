#!/bin/bash
## Configuration values for SLURM job submission.
## One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --job-name=getVariantCounts        # friendly name for job.
#SBATCH --nodes=1                      # ensure cpus are on one node
#SBATCH --ntasks=1                     # run a single task
#SBATCH --cpus-per-task=4              # number of cpus/threads requested.
#SBATCH --mem=16gb                     # memory requested.
#SBATCH --partition=20                 # partition (queue) to use
#SBATCH --output=/lab/jain_imaging/Kelsey/Slurm_out/getVariantCounts_%j.out             # name of output file.  %j is jobid
#SBATCH --error=/lab/jain_imaging/Kelsey/Slurm_out/getVariantCounts_%j.err
#SBATCH --mail-type=END                # send email on job start/finish.
#SBATCH --mail-user=farenhem@wi.mit.edu

path=/lab/jain_imaging/Kelsey/Sequencing/20210903_NovaSeq/ANALYSIS
pileuppath=${path}/variant_pileups
codepath=${path}/code
sample=$1

for chrom in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
do
for strand in pos neg 
do 

gunzip ${sample}_${chrom}_vars.pileup.gz
bedmap --echo --echo-map-id-uniq ${sample}_${chrom}_vars.pileup ${refpath}/hg38knownCanonicalGenes_StrandOnly.bed | sed 's/|/\t/g' - > ${sample}_${chrom}_vars_strand.pileup
mv ${sample}_${chrom}_vars_strand.pileup ${sample}_${chrom}_vars.pileup

done
done

python -u ${codepath}/compile_variant_counts.py --path ${pileuppath} --sample ${sample}