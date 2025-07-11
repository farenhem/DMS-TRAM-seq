#!/bin/sh
#SBATCH --job-name=join_windows
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --partition=20

# Adds header and combine chromosomes in specified order, both for tsv and csv files. Then, apply p-val correction. 

windowSize=$1

cd ../windows

echo "gene,ID,type,chr,strand,region,exon,start,end,x,delta,p,r_con1,r_con2,r_between,con1_avg,con2_avg,avg_pval,con1_gini,con2_gini,gini_pval,n,len" | cat - ${windowSize}nt_windows_step_chr1.csv > ${windowSize}nt_windows_step_temp1.csv
echo "gene,ID,type,chr,strand,region,exon,start,end,x,delta,p,r_con1,r_con2,r_between,con1_avg,con2_avg,avg_pval,con1_gini,con2_gini,gini_pval,n,len" | cat - ${windowSize}nt_windows_slide_chr1.csv > ${windowSize}nt_windows_slide_temp1.csv
echo "gene\tID\ttype\tchr\tstrand\tregion\texon\tstart\tend\tx\tdelta\tp\tr_con1\tr_con2\tr_between\tcon1_avg\tcon2_avg\tavg_pval\tcon1_gini\tcon2_gini\tgini_pval\tn\tlen" | cat - ${windowSize}nt_windows_step_chr1.txt > ${windowSize}nt_windows_step_temp1.txt
echo "gene\tID\ttype\tchr\tstrand\tregion\texon\tstart\tend\tx\tdelta\tp\tr_con1\tr_con2\tr_between\tcon1_avg\tcon2_avg\tavg_pval\tcon1_gini\tcon2_gini\tgini_pval\tn\tlen" | cat - ${windowSize}nt_windows_slide_chr1.txt > ${windowSize}nt_windows_slide_temp1.txt


for chrom in chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
do 

cat ${windowSize}nt_windows_step_temp1.csv ${windowSize}nt_windows_step_${chrom}.csv > ${windowSize}nt_windows_step_temp2.csv
mv ${windowSize}nt_windows_step_temp2.csv ${windowSize}nt_windows_step_temp1.csv
cat ${windowSize}nt_windows_slide_temp1.csv ${windowSize}nt_windows_slide_${chrom}.csv > ${windowSize}nt_windows_slide_temp2.csv
mv ${windowSize}nt_windows_slide_temp2.csv ${windowSize}nt_windows_slide_temp1.csv

cat ${windowSize}nt_windows_step_temp1.txt ${windowSize}nt_windows_step_${chrom}.txt > ${windowSize}nt_windows_step_temp2.txt
mv ${windowSize}nt_windows_step_temp2.txt ${windowSize}nt_windows_step_temp1.txt
cat ${windowSize}nt_windows_slide_temp1.txt ${windowSize}nt_windows_slide_${chrom}.txt > ${windowSize}nt_windows_slide_temp2.txt
mv ${windowSize}nt_windows_slide_temp2.txt ${windowSize}nt_windows_slide_temp1.txt

done

mv ${windowSize}nt_windows_step_temp1.csv ${windowSize}nt_windows_step.csv
mv ${windowSize}nt_windows_slide_temp1.csv ${windowSize}nt_windows_slide.csv
mv ${windowSize}nt_windows_step_temp1.txt ${windowSize}nt_windows_step.txt
mv ${windowSize}nt_windows_slide_temp1.txt ${windowSize}nt_windows_slide.txt

mkdir -p chroms
mv ${windowSize}nt_*_chr*.csv chroms/
mv ${windowSize}nt_*_chr*.txt chroms/

python windows_padj_BH.py ${windowSize}nt_windows_step.csv &
python windows_padj_BH.py ${windowSize}nt_windows_slide.csv &

wait

mv ${windowSize}nt_windows_step_padj.csv ${windowSize}nt_windows_step.csv
mv ${windowSize}nt_windows_slide_padj.csv ${windowSize}nt_windows_slide.csv
mv ${windowSize}nt_windows_step_padj.txt ${windowSize}nt_windows_step.txt
mv ${windowSize}nt_windows_slide_padj.txt ${windowSize}nt_windows_slide.txt
