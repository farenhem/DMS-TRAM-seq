#!/bin/sh
#SBATCH --job-name=bedmap_CLIP
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --partition=20
#SBATCH --output /lab/jain_imaging/Kelsey/Slurm_out/window_CLIP_%j.out
#SBATCH --error /lab/jain_imaging/Kelsey/Slurm_out/window_CLIP_%j.err
#Add header and combine chromosomes in specified order, both for tsv and csv files. Then, apply p-val correction. 

windowSize=$1

path=/lab/jain_imaging/Kelsey/Sequencing/20210903_NovaSeq/ANALYSIS
refpath=${path}/reference_annot

awk -v FS='\t' -v OFS='\t' '{print $4, $8, $9, $2, $1, $5, $3, $6, $7, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24}' ${windowSize}nt_windows_step.txt | sed '1d' > ${windowSize}nt_windows_step_temp.bed
awk -v FS='\t' -v OFS='\t' '{print $4, $8, $9, $2, $1, $5, $3, $6, $7, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24}' ${windowSize}nt_windows_slide.txt | sed '1d' > ${windowSize}nt_windows_slide_temp.bed



bedtools intersect -u -a ${refpath}/eCLIP_HepG2_abbr.bed -b ${refpath}/canon_GRCh38.106_all_exons.bed > ${refpath}/HepG2_filtered.bed
bedtools intersect -u -a ${refpath}/eCLIP_K562_abbr.bed -b ${refpath}/canon_GRCh38.106_all_exons.bed > ${refpath}/K562_filtered.bed


sort-bed ${windowSize}nt_windows_step_temp.bed > ${windowSize}nt_windows_step_sorted.bed
sort-bed ${windowSize}nt_windows_slide_temp.bed > ${windowSize}nt_windows_slide_sorted.bed

#annotate with CLIP peaks overlapping with the window
bedmap --range 0 --echo --echo-map-id ${windowSize}nt_windows_step_sorted.bed ${refpath}/HepG2_filtered.bed | sed 's/|/\t/g' - > ${windowSize}nt_windows_step_sorted_HepG2.bed
bedmap --range 0 --echo --echo-map-id ${windowSize}nt_windows_step_sorted_HepG2.bed ${refpath}/K562_filtered.bed | sed 's/|/\t/g' - > ${windowSize}nt_windows_step_CLIP.bed
bedmap --range 0 --echo --echo-map-id ${windowSize}nt_windows_slide_sorted.bed ${refpath}/HepG2_filtered.bed | sed 's/|/\t/g' - > ${windowSize}nt_windows_slide_sorted_HepG2.bed
bedmap --range 0 --echo --echo-map-id ${windowSize}nt_windows_slide_sorted_HepG2.bed ${refpath}/K562_filtered.bed | sed 's/|/\t/g' - > ${windowSize}nt_windows_slide_CLIP.bed

rm ${windowSize}nt_windows_step_temp.bed ${windowSize}nt_windows_step_sorted.bed ${windowSize}nt_windows_step_sorted_HepG2.bed
rm ${windowSize}nt_windows_slide_temp.bed ${windowSize}nt_windows_slide_sorted.bed ${windowSize}nt_windows_slide_sorted_HepG2.bed

sed 's/,/;/g' ${windowSize}nt_windows_step_CLIP.bed | awk 'BEGIN { FS="\t"; OFS="," } {$1=$1; print}' - > ${windowSize}nt_windows_step_CLIP_noheader.csv
sed 's/,/;/g' ${windowSize}nt_windows_slide_CLIP.bed | awk 'BEGIN { FS="\t"; OFS="," } {$1=$1; print}' - > ${windowSize}nt_windows_slide_CLIP_noheader.csv

echo "chr,start,end,ID,gene,strand,type,region,exon,x,delta,p,padj,r_con1,r_con2,r_between,con1_avg,con2_avg,avg_pval,con1_gini,con2_gini,gini_pval,n,len,HepG2,K562" | cat - ${windowSize}nt_windows_step_CLIP_noheader.csv > ${windowSize}nt_windows_step_CLIP.csv
echo "chr,start,end,ID,gene,strand,type,region,exon,x,delta,p,padj,r_con1,r_con2,r_between,con1_avg,con2_avg,avg_pval,con1_gini,con2_gini,gini_pval,n,len,HepG2,K562" | cat - ${windowSize}nt_windows_slide_CLIP_noheader.csv > ${windowSize}nt_windows_slide_CLIP.csv

rm ${windowSize}nt_windows_step_CLIP_noheader.csv ${windowSize}nt_windows_slide_CLIP_noheader.csv

