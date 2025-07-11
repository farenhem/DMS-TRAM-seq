#!/bin/bash

# This script is a control script for genome-wide analysis of DMS-MaPseq
# data in different conditions, starting with BAM files (i.e. mapped reads)
# as input.

# This is step 3 of 8, where DMS mismatch rates are normalized using control
# (i.e. non DMS-treated) rates following the subtraction convention.

# Recall that the control samples span a genomic footprint where they have
# coverage above a specified threshold, implying the same footprint for the
# normalized rates that are calculated here.

jobpath=../jobs
pileuppath=../pileups
bedGraphpath=../bedGraph

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

# Pair up the DMS experiments with non-DMS controls.
declare -A pairs=( [${DMS1}]=${ctl1} [${DMS2}]=${ctl2} [${DMS3}]=${ctl3} [${sDMS1}]=${sctl1} [${sDMS2}]=${sctl2} [${sDMS3}]=${sctl3} )

# Normalize using the convention of subtracting the control mismatch rate from
# the DMS mismatch rate.
for treat in "${!pairs[@]}"
do
    ctl=${pairs[${treat}]}
    for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
    do
	
	for strand in pos neg
	do

cat - << _Eod1_ > ${jobpath}/intCoverage_${treat}_${chr}_${strand}.slurm
#!/bin/bash
#SBATCH --job-name=int_cov_${treat}_${chr}_${strand}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --partition=20
#SBATCH --output ${jobpath}/int_cov_${treat}_${chr}_${strand}.out

cd ${bedGraphpath}

bedtools intersect -wb -wa -a ${treat}_coverage_${chr}_${strand}.bedGraph -b ${ctl}_coverage_${chr}_${strand}.bedGraph | awk '{OFS="\t"}{print \$1,\$2,\$3,\$4}' > ${treat}_coverage_norm_${chr}_${strand}.bedGraph
rm ${treat}_coverage_${chr}_${strand}.bedGraph ${ctl}_coverage_${chr}_${strand}.bedGraph
_Eod1_
	    
cat - << _Eod1_ > ${jobpath}/normRate_${treat}_${chr}_${strand}.slurm
#!/bin/bash
#SBATCH --job-name=norm_${treat}_${chr}_${strand}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --partition=20
#SBATCH --output ${jobpath}/norm_${treat}_${chr}_${strand}.out

cd ${bedGraphpath}

bedtools intersect -wb -wa -a ${treat}_mmRate_${chr}_${strand}.bedGraph -b ${ctl}_mmRate_${chr}_${strand}.bedGraph | awk '{max=0}{if (\$4>\$8) max=\$4-\$8}{OFS="\t"}{print \$1,\$2,\$3,max}' > ${treat}_mmRate_norm_${chr}_${strand}.bedGraph
rm ${treat}_mmRate_${chr}_${strand}.bedGraph ${ctl}_mmRate_${chr}_${strand}.bedGraph
_Eod1_

	done
    done
done

for treat in ${DMS1} ${DMS2} ${DMS3} ${sDMS1} ${sDMS2} ${sDMS3}
do
    for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
    do
	for strand in pos neg
	do
	    sbatch ${jobpath}/intCoverage_${treat}_${chr}_${strand}.slurm
	    sbatch ${jobpath}/normRate_${treat}_${chr}_${strand}.slurm
	done
    done
done

exit
