#!/bin/sh
#SBATCH --job-name=GWSA4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --partition=20

# This script is a control script for genome-wide analysis of DMS-MaPseq
# data in different conditions, starting with BAM files (i.e. mapped reads)
# as input.

# This is step 4 of 8, where the files for per-chromosome normalized mismatch
# rates are combined into single bedGraph files for each condition.

# List experiments from non-stress conditions.
ctl1=C1
ctl2=C2
ctl3=C3
DMS1=D1
DMS2=D2
DMS3=D3

# List experiments from stress conditions.
sctl1=A1
sctl2=A2
sctl3=A3
sDMS1=AD1
sDMS2=AD2
sDMS3=AD3

# After the coverage and mismatch rates have been tabulated and normalized
# for each chromosome, reassemple the genome.

jobpath=../jobs
pileuppath=../pileups
bedGraphpath=../bedGraph
cd ${bedGraphpath}

for experiment in ${DMS1} ${DMS2} ${DMS3} ${sDMS1} ${sDMS2} ${sDMS3}
do
    for assay in coverage mmRate
    do		 
	for strand in pos neg
	do
	    cat ${experiment}_${assay}_norm_chr1_${strand}.bedGraph ${experiment}_${assay}_norm_chr2_${strand}.bedGraph > temp_si.bedGraph
	    rm ${experiment}_${assay}_norm_chr1_${strand}.bedGraph ${experiment}_${assay}_norm_chr2_${strand}.bedGraph
            for chr in chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
            do
		cat temp_si.bedGraph ${experiment}_${assay}_norm_${chr}_${strand}.bedGraph > temp_sj.bedGraph
		mv temp_sj.bedGraph temp_si.bedGraph
		rm ${experiment}_${assay}_norm_${chr}_${strand}.bedGraph
	    done
	    mv temp_si.bedGraph ${experiment}_${assay}_${strand}.bedGraph
	done
    done
done

exit
