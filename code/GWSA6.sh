#!/bin/sh

# This script is a control script for genome-wide analysis of DMS-MaPseq
# data in different conditions, starting with BAM files (i.e. mapped reads)
# as input.

# This is step 6 of 8, where sites that were retained in step 5 are used
# to set common sites for all samples.  After this step, sites are tabulated
# (on positive and negative strands, respectively) are common to all samples.

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

for treat in ${DMS1} ${DMS2} ${DMS3} ${sDMS1} ${sDMS2} ${sDMS3}
do
    for assay in coverage mmRate
    do
	for strand in pos neg
	do
	    
cat - << _Eod1_ > ${jobpath}/setCommon_${treat}_${assay}_${strand}.slurm
#!/bin/bash
#SBATCH --job-name=setCommon_${treat}_${assay}_${strand}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --partition=20
#SBATCH --output ${jobpath}/setCommon_${treat}_${assay}_${strand}.out

cd ${bedGraphpath}

mv ${treat}_${assay}_${strand}.bedGraph temp_${treat}_${assay}_${strand}.bedGraph
bedtools intersect -wb -wa -a temp_${treat}_${assay}_${strand}.bedGraph -b commonSites_coverage_${strand}.bedGraph | awk '{OFS="\t"}{print \$1,\$2,\$3,\$4}' > ${treat}_${assay}_${strand}.bedGraph
rm temp_${treat}_${assay}_${strand}.bedGraph

_Eod1_

	done
    done
done

for treat in ${DMS1} ${DMS2} ${DMS3} ${sDMS1} ${sDMS2} ${sDMS3}
do
    for assay in coverage mmRate
    do
	for strand in pos neg
	do
	    sbatch ${jobpath}/setCommon_${treat}_${assay}_${strand}.slurm
	    #sleep 1
	done
    done
done

exit
