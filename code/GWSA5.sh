#!/bin/sh
#SBATCH --job-name=GWSA5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --partition=20
#SBATCH --output /lab/jain_imaging/Kelsey/Sequencing/20210903_NovaSeq/ANALYSIS/jobs/GWSA5_%j.out

# This script is a control script for genome-wide analysis of DMS-MaPseq
# data in different conditions, starting with BAM files (i.e. mapped reads)
# as input.

# This is step 5 of 8, where sites are retained only if they exist (i.e.
# have already passed the coverage threshold condition) in ALL samples.
# Note that the earlier normalization step has already partly imposed this
# condition, as retained (i.e. normalized) sites must have passed the
# threshold criterion in each control/DMS-treated sample pair.

path=/lab/jain_imaging/Kelsey/Sequencing/20210903_NovaSeq
jobpath=${path}/ANALYSIS/jobs
bedGraphpath=${path}/ANALYSIS/bedGraph

# List experiments from non-stress conditions.
DMS1=D1
DMS2=D2
DMS3=D3

# List experiments from stress conditions.
sDMS1=AD1
sDMS2=AD2
sDMS3=AD3

cd ${bedGraphpath}

for strand in pos neg
do	      
    bedtools intersect -wb -wa -a ${DMS1}_coverage_${strand}.bedGraph -b ${DMS2}_coverage_${strand}.bedGraph | awk '{OFS="\t"}{print $1,$2,$3,$4}' > temp.bedGraph
    for experiment in ${DMS3} ${sDMS1} ${sDMS2} ${sDMS3}
    do
	bedtools intersect -wb -wa -a temp.bedGraph -b ${experiment}_coverage_${strand}.bedGraph | awk '{OFS="\t"}{print $1,$2,$3,$4}' > temp2.bedGraph
	mv temp2.bedGraph temp.bedGraph
    done
    mv temp.bedGraph commonSites_coverage_${strand}.bedGraph
done

exit
