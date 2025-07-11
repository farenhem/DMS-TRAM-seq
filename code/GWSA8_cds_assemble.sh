#!/bin/sh
#SBATCH --job-name=GWSA8_cds
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --partition=20

# This script is a control script for genome-wide analysis of DMS-MaPseq
# data in different conditions, starting with BAM files (i.e. mapped reads)
# as input.

# This is step 8 of 8, where per-base normalized, Winsorized mismatch rates
# from DMS-treated samples are used to compute per-element (e.g. 3' UTRs,
# exons, introns, canonical transcripts, etc.) averages for coverage
# mismatch rates, entropy and Gini index.

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

codepath=../code
annotpath=../reference_annot
bedGraphpath=../bedGraph
cd ${bedGraphpath}

# For a set of reference annotations (e.g. complete transcripts, individual
# exons, 3' UTRs, introns etc.), collect per-nucleotide coverage and
# mismatch rates.

# Set the strand-specific references.
reference=canon_GRCh38.106_cds

# First list all of the candidate sites in each region.
for treat in ${DMS1} ${DMS2} ${DMS3} ${sDMS1} ${sDMS2} ${sDMS3}
do
    for assay in coverage mmRate
    do
	for strand in pos neg
	do
	    bedtools intersect -wb -wa -a ${annotpath}/${reference}_${strand}.bed -b ${treat}_${assay}_${strand}.bedGraph > ${treat}_${assay}_${strand}_cds.bed
	    bedtools groupby -i ${treat}_${assay}_${strand}_cds.bed -g 1,2,3,4,5,6 -c 10 -o collapse > ${treat}_${assay}_${strand}_cds_collapse.txt
	    awk '{OFS="\t"}{print $1,$2,$3,$4,$7,$6}' ${treat}_${assay}_${strand}_cds_collapse.txt > ${treat}_${assay}_${strand}_cds_collapse.bed	    
	    rm ${treat}_${assay}_${strand}_cds.bed ${treat}_${assay}_${strand}_cds_collapse.txt
	done
	cat ${treat}_${assay}_pos_cds_collapse.bed ${treat}_${assay}_neg_cds_collapse.bed > ${treat}_${assay}_cds_collapsed.bed
	sort -k 1,1 -k4,4 -k2,2n ${treat}_${assay}_cds_collapsed.bed > ${treat}_${assay}_cds_collapse.bed
	rm ${treat}_${assay}_pos_cds_collapse.bed ${treat}_${assay}_neg_cds_collapse.bed ${treat}_${assay}_cds_collapsed.bed
    done
done

# Next, aggregate regions into transcripts. Note, it should be possible
# to combine the above and below loops into a single step.

for treat in ${DMS1} ${DMS2} ${DMS3} ${sDMS1} ${sDMS2} ${sDMS3}
do
    for assay in coverage mmRate
    do
	bedtools groupby -i ${treat}_${assay}_cds_collapse.bed -g 4 -c 1 -o first | awk '{print $2}' > ${treat}_${assay}_cds_collapseGenes_chr.tmp
	bedtools groupby -i ${treat}_${assay}_cds_collapse.bed -g 4 -c 2 -o min | awk '{print $2}' > ${treat}_${assay}_cds_collapseGenes_start.tmp
	bedtools groupby -i ${treat}_${assay}_cds_collapse.bed -g 4 -c 3 -o max | awk '{print $2}' > ${treat}_${assay}_cds_collapseGenes_end.tmp	
	bedtools groupby -i ${treat}_${assay}_cds_collapse.bed -g 4 -c 6 -o first | awk '{print $2}' > ${treat}_${assay}_cds_collapseGenes_strand.tmp
	bedtools groupby -i ${treat}_${assay}_cds_collapse.bed -g 4 -c 5 -o collapse > ${treat}_${assay}_cds_collapseGenes_sig.tmp
        pr -m -t -J ${treat}_${assay}_cds_collapseGenes_chr.tmp ${treat}_${assay}_cds_collapseGenes_start.tmp ${treat}_${assay}_cds_collapseGenes_end.tmp ${treat}_${assay}_cds_collapseGenes_sig.tmp ${treat}_${assay}_cds_collapseGenes_strand.tmp | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6}' > ${treat}_${assay}_cds_collapseGenes.tmp
	sort -k 1,1 -k 2,2n ${treat}_${assay}_cds_collapseGenes.tmp > ${treat}_${assay}_cds_collapseGenes.bed
	rm ${treat}_${assay}_cds_collapseGenes*.tmp
    done
done

# For any canonical transcripts that have been shortened here due to exons not
# being expressed up to the specified threshold(s), update the coordinates to
# reflect the full length transcript.

for treat in ${DMS1} ${DMS2} ${DMS3} ${sDMS1} ${sDMS2} ${sDMS3}
do
    for assay in coverage mmRate
    do
	mv ${treat}_${assay}_cds_collapseGenes.bed ${treat}_${assay}_cds_collapseGenes.tmp
	python ${codepath}/updateBed.py --ref ${annotpath}/canon_GRCh38.106_cds_transcripts.bed --query ${treat}_${assay}_cds_collapseGenes.tmp --out ${treat}_${assay}_cds_collapseGenes.bed
	rm ${treat}_${assay}_cds_collapseGenes.tmp
    done
done

for assay in coverage mmRate
do
    awk '{OFS="\t"}{print $1,$2,$3,$4,$6,$5}' ${DMS1}_${assay}_cds_collapseGenes.bed > ${DMS1}_cds.tmp
    awk '{print $5}' ${DMS2}_${assay}_cds_collapseGenes.bed > ${DMS2}_cds.tmp
    awk '{print $5}' ${DMS3}_${assay}_cds_collapseGenes.bed > ${DMS3}_cds.tmp
    awk '{print $5}' ${sDMS1}_${assay}_cds_collapseGenes.bed > ${sDMS1}_cds.tmp
    awk '{print $5}' ${sDMS2}_${assay}_cds_collapseGenes.bed > ${sDMS2}_cds.tmp
    awk '{print $5}' ${sDMS3}_${assay}_cds_collapseGenes.bed > ${sDMS3}_cds.tmp
    pr -m -t -J ${DMS1}_cds.tmp ${DMS2}_cds.tmp ${DMS3}_cds.tmp ${sDMS1}_cds.tmp ${sDMS2}_cds.tmp ${sDMS3}_cds.tmp | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' > ${assay}CDSAllGenes.tmp
    cat ${codepath}/${assay}Head6.txt ${assay}CDSAllGenes.tmp > ${assay}CDSAllGenes.txt
    rm -f ${DMS1}_cds.tmp ${DMS2}_cds.tmp ${DMS3}_cds.tmp ${sDMS1}_cds.tmp ${sDMS2}_cds.tmp ${sDMS3}_cds.tmp
    rm -f ${DMS1}_${assay}_cds_collapse*.bed ${DMS2}_${assay}_cds_collapse*.bed ${DMS3}_${assay}_cds_collapse*.bed ${sDMS1}_${assay}_cds_collapse*.bed ${sDMS2}_${assay}_cds_collapse*.bed ${sDMS3}_${assay}_cds_collapse*.bed
    rm -f ${assay}CDSAllGenes.tmp
done

python ${codepath}/process_replicate_signal.py --covin coverageCDSAllGenes.txt --mmin mmRateCDSAllGenes.txt --outfile processedExperimentsCDS.csv

exit
