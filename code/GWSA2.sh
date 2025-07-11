#!/bin/sh

# This script is a control script for genome-wide analysis of DMS-MaPseq
# data in different conditions, starting with BAM files (i.e. mapped reads)
# as input.

# This is step 2 of 8, where a threshold constraint on coverage is applied
# genome-wide to all samples.  At the end of this step, the genomic footprint
# for the control samples is comprised of only bases where coverage is above
# the specified threshold.

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

# Make sure that coverage is above a defined threshold in the control
# and treatment samples.
threshold=100
jobpath=../jobs
pileuppath=../pileups
bedGraphpath=../bedGraph

mkdir -p ${bedGraphpath}

for experiment in ${ctl1} ${ctl2} ${ctl3} ${DMS1} ${DMS2} ${DMS3} ${sctl1} ${sctl2} ${sctl3} ${sDMS1} ${sDMS2} ${sDMS3}
do
    for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
    do

cat - << _Eod1_ > ${jobpath}/filtCoverage_${experiment}_${chr}.slurm
#!/bin/bash
#SBATCH --job-name=filter_cov_${experiment}_${chr}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --partition=20
#SBATCH --output ${jobpath}/filter_cov_${experiment}_${chr}.out

cd ${pileuppath}

awk -v cut=${threshold} '{ if (\$5+\$6 >= cut) {{OFS="\t"}{print \$1,\$2-1,\$2,\$5+\$6}}}' ${experiment}_${chr}_pos.pileup | grep -v "pos" > ${bedGraphpath}/${experiment}_coverage_${chr}_pos.bedGraph
awk -v cut=${threshold} '{ if (\$5+\$6 >= cut) {{OFS="\t"}{print \$1,\$2-1,\$2,\$5+\$6}}}' ${experiment}_${chr}_neg.pileup | grep -v "pos" > ${bedGraphpath}/${experiment}_coverage_${chr}_neg.bedGraph
awk -v cut=${threshold} '{ if (\$5+\$6 >= cut) {{OFS="\t"}{print \$1,\$2-1,\$2,\$7}}}' ${experiment}_${chr}_pos.pileup | grep -v "pos" > ${bedGraphpath}/${experiment}_mmRate_${chr}_pos.bedGraph
awk -v cut=${threshold} '{ if (\$5+\$6 >= cut) {{OFS="\t"}{print \$1,\$2-1,\$2,\$7}}}' ${experiment}_${chr}_neg.pileup | grep -v "pos" > ${bedGraphpath}/${experiment}_mmRate_${chr}_neg.bedGraph
gzip ${experiment}_${chr}_pos.pileup
gzip ${experiment}_${chr}_neg.pileup
_Eod1_
	
    done
done

for experiment in ${ctl1} ${ctl2} ${ctl3} ${DMS1} ${DMS2} ${DMS3} ${sctl1} ${sctl2} ${sctl3} ${sDMS1} ${sDMS2} ${sDMS3}
do
    for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
    do
	sbatch ${jobpath}/filtCoverage_${experiment}_${chr}.slurm
    done
done

exit
