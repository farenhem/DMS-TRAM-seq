#!/bin/bash
## Configuration values for SLURM job submission.
## One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --job-name=window_gen        # friendly name for job.
#SBATCH --nodes=1                      # ensure cpus are on one node
#SBATCH --ntasks=25                    # run 25 tasks (one per chrom)
#SBATCH --cpus-per-task=1              # number of cpus/threads requested.
#SBATCH --mem=64gb                     # memory requested.
#SBATCH --partition=20                 # partition (queue) to use
#SBATCH --output=/lab/jain_imaging/Kelsey/Slurm_out/window_gen_%j.out             # name of output file.  %j is jobid
#SBATCH --error=/lab/jain_imaging/Kelsey/Slurm_out/window_gen_%j.err
#SBATCH --mail-type=END                # send email on job start/finish.
#SBATCH --mail-user=farenhem@wi.mit.edu

path=/lab/jain_imaging/Kelsey/Sequencing/20210903_NovaSeq/ANALYSIS
bedGraphpath=${path}/bedGraph
codepath=${path}/code
refpath=${path}/reference_annot
winpath=${path}/set_windows
winsize=$1

cd ${bedGraphpath}


<<comment
###if necessary, separate bedGraph files by chromosome
for sample in D1 D2 D3 AD1 AD2 AD3
do
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY 
do
for strand in pos neg
do
awk -v chrom="$chr" '$1==chrom' ${sample}_mmRate_${strand}.bedGraph > chrom/${sample}_${chr}_mmRate_${strand}.bedGraph
awk -v chrom="$chr" '$1==chrom' ${sample}_coverage_${strand}.bedGraph > chrom/${sample}_${chr}_coverage_${strand}.bedGraph
done
done
done




###generate sliding and step windows; for sliding windows, also removes overlapping windows, keeping those with higher delta values
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX #chrY 
do
srun --exclusive --mem 8G -n1 python -u ${codepath}/window_generate_1.py --path ${path} --chr ${chr} --size ${winsize} --condition1 D --condition2 AD & #ampersand parallelizes tasks
done

wait #all chromosomes run
comment
#remove CLIP peaks where there is no coverage; this prevents including peaks where there are is no data later on

for strand in pos neg 
do
bedtools intersect -u -a ${refpath}/eCLIP_HepG2_abbr.bed -b ${bedGraphpath}/D1_coverage_${strand}.bedGraph > ${refpath}/HepG2_filtered_${strand}.bed
bedtools intersect -u -a ${refpath}/eCLIP_K562_abbr.bed -b ${bedGraphpath}/D1_coverage_${strand}.bedGraph > ${refpath}/K562_filtered_${strand}.bed
done

###add strand-specific transcript annotation, UTR/exon info, and eCLIP peaks from ENCODE datasets (HepG2 and K562), which include many RBPs
cd ${winpath}
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX #chrY  
do
for type in SLIDE STEP
do
for strand in pos neg
do
for comparison in D_AD 
do
sort-bed chrom/${comparison}_${chr}_${winsize}bp_deltas_${type}_${strand}.txt > chrom/${comparison}_${chr}_${winsize}bp_deltas_${type}_sorted_${strand}.txt

#annotate with transcript info, as well as information on the specific region (ie CDS vs UTR)
bedmap --delim '\t' --echo --echo-map-id-uniq chrom/${comparison}_${chr}_${winsize}bp_deltas_${type}_sorted_${strand}.txt ${refpath}/GENCODEv43_knownCanonical_allRegions.bed > chrom/${comparison}_${chr}_${winsize}bp_deltas_${type}_annot_${strand}_IDonly.txt
bedmap --delim '\t' --echo --echo-map-id chrom/${comparison}_${chr}_${winsize}bp_deltas_${type}_annot_${strand}_IDonly.txt ${refpath}/GENCODEv43_knownCanonical_allRegions_regiononly.bed > chrom/${comparison}_${chr}_${winsize}bp_deltas_${type}_annot_${strand}_temp.txt

#annotate with CLIP peaks overlapping with the window
bedmap --delim '\t' --range 0 --echo --echo-map-id chrom/${comparison}_${chr}_${winsize}bp_deltas_${type}_annot_${strand}_temp.txt ${refpath}/HepG2_filtered_${strand}.bed > chrom/${comparison}_${chr}_${winsize}bp_deltas_${type}_annot_${strand}_temp2.txt
bedmap --delim '\t' --range 0 --echo --echo-map-id chrom/${comparison}_${chr}_${winsize}bp_deltas_${type}_annot_${strand}_temp2.txt ${refpath}/K562_filtered_${strand}.bed > chrom/${comparison}_${chr}_${winsize}bp_deltas_${type}_annot_${strand}.txt

#rm chrom/${comparison}_${chr}_${winsize}bp_deltas_${type}_${strand}.txt
rm chrom/${comparison}_${chr}_${winsize}bp_deltas_${type}_sorted_${strand}.txt  
rm chrom/${comparison}_${chr}_${winsize}bp_deltas_${type}_annot_${strand}_IDonly.txt
rm chrom/${comparison}_${chr}_${winsize}bp_deltas_${type}_annot_${strand}_temp.txt
rm chrom/${comparison}_${chr}_${winsize}bp_deltas_${type}_annot_${strand}_temp2.txt
done #comparison loop
done #strand loop
done #type loop
done #chr loop






###add info based on transcript ID, such as transcript type and stress granule localization (from Parker paper)
for comparison in D_AD 
do

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX #chrY  
do
srun --exclusive --mem 3G -n1 python -u ${codepath}/window_generate_2.py --path ${path} --chr ${chr} --size ${winsize} --comparison ${comparison} & 
done 


wait 

rm chrom/${comparison}_${chr}_${winsize}bp_deltas_*_annot_*.txt

cd ${winpath}

cat chrom/${comparison}_chr*_${winsize}bp_SLIDE_deltas.csv | sort -gr -k8 -t, > ${comparison}_${winsize}bp_SLIDE_deltas_noheader.csv
cat chrom/${comparison}_chr*_${winsize}bp_STEP_deltas.csv | sort -gr -k8 -t, > ${comparison}_${winsize}bp_STEP_deltas_noheader.csv

echo -e "chr,start,end,ID,gene,type,exon,strand,delta,p,con1_avg,con2_avg,avg_pval,con1_gini,con2_gini,gini_pval,HepG2,K562,SG_FC,SG_localization" | cat - ${comparison}_${winsize}bp_SLIDE_deltas_noheader.csv > ${comparison}_${winsize}bp_SLIDE_deltas.csv #pos and neg already joined in python step 
echo -e "chr,start,end,ID,gene,type,exon,strand,delta,p,con1_avg,con2_avg,avg_pval,con1_gini,con2_gini,gini_pval,HepG2,K562,SG_FC,SG_localization" | cat - ${comparison}_${winsize}bp_STEP_deltas_noheader.csv > ${comparison}_${winsize}bp_STEP_deltas.csv #pos and neg already joined in python step

rm ${comparison}_${winsize}bp_SLIDE_deltas_noheader.csv
rm chrom/${comparison}_chr*_${winsize}bp_*_deltas.csv 
done



###once the compiled (all chromosome) csv is made, compile a list of transcripts containing the most significantly-changed windows. Also applies p-val correction and rewrites csvs above.
for comparison in D_AD 
do
srun --exclusive --mem 3G -n1 python -u ${codepath}/window_generate_3.py --path ${path} --size ${winsize} --comparison ${comparison} & 
done

wait
