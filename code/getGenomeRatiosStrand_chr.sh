#!/bin/sh

refpath=../reference_annot/fasta_hg38
datapath=../BAM_primary
jobpath=../jobs
codepath=../code
pileuppath=../pileups
sizes=../reference_annot/fasta_hg38/chromInfo.txt
bamfileroot=$1

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
do

cat - << _Eod1_ > ${jobpath}/getRatiosStrand_${bamfileroot}_${chr}.slurm
#!/bin/bash
#SBATCH --job-name=ratios_${bamfileroot}_${chr}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --partition=20
#SBATCH --output ${jobpath}/ratios_${bamfileroot}_${chr}.out

cd ${pileuppath}

# Extract chromosome-specific bam file.
samtools view ${datapath}/${bamfileroot}_primary.bam ${chr} -b > ${datapath}/${bamfileroot}_${chr}.bam

# Index chromosome-specific bam file.
samtools index ${datapath}/${bamfileroot}_${chr}.bam ${datapath}/${bamfileroot}_${chr}.bai

# Get raw pileup file from bcftools, cleaning up the working dir. along the way.
bcftools mpileup -d 1000000 -q 20 -Q 20 -a format/ad,info/ad -Ov -f ${refpath}/${chr}.fa ${datapath}/${bamfileroot}_${chr}.bam > ${bamfileroot}_${chr}.vcf
rm ${datapath}/${bamfileroot}_${chr}.bam ${datapath}/${bamfileroot}_${chr}.bai

bcftools norm -f ${refpath}/${chr}.fa ${bamfileroot}_${chr}.vcf > ${bamfileroot}_${chr}_norm.vcf
rm -f ${bamfileroot}_${chr}.vcf

bcftools query -f '%CHROM  %POS  %REF  %TYPE  %ALT{0}  %IMF  %DP  %AD\n' ${bamfileroot}_${chr}_norm.vcf | awk '{OFS="\t"}{print \$1"_"\$2,\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8}' | sed 's/<\*>/R/g' > ${bamfileroot}_${chr}_norm.pileupt
rm -f ${bamfileroot}_${chr}_norm.vcf

# Process raw pileup file to get ratios of interest.
python ${codepath}/process_pileup_strand_select.py --pileup ${bamfileroot}_${chr}_norm.pileupt --outfile1 ${bamfileroot}_${chr}_pos.pileup --outfile2 ${bamfileroot}_${chr}_neg.pileup --outfile3 ${bamfileroot}_${chr}_vars.pileup
rm -f ${bamfileroot}_${chr}_norm.pileupt
gzip ${bamfileroot}_${chr}_vars.pileup

_Eod1_

done

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
do
    sbatch ${jobpath}/getRatiosStrand_${bamfileroot}_${chr}.slurm
    #sleep 1
done

exit
