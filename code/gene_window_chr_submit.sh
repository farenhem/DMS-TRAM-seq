#!/bin/sh
##SBATCH --job-name=gene_windows
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --partition=20

chrom=$1
winSize=$2
bgpath=../bedGraph/chrom
windowpath=../gene_windows
ref=../reference_annot/canon_GRCh38.106_all_exons.bed


sample1=D1
sample2=D2
sample3=D3
sample4=AD1
sample5=AD2
sample6=AD3



cd ${windowpath}
mkdir -p gene_annot
mkdir -p gene_bedGraph

python -u gene_window_gen.py --chr ${chrom} --bgpath ${bgpath} --out ${windowpath} --ref ${ref} --winsize ${winSize}
