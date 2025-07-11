# DMS-TRAM-seq
This software will fully analyze DMS-TRAM-seq data, from raw FASTQ files through the analyzed mutational profiles. Please note this work is under review and this pipeline will be updated. See our preprint here: https://www.biorxiv.org/content/10.1101/2025.06.05.658101v1

# Installation

Please note, as written, this code is designed to function on a SLURM computing cluster. Modification will be needed to run on another cluster or on a single computer. These updates are in progress.

Installing this software and the required reference data should take about 20 minutes. 

## Pre-requisites

All software needed for this pipeline is described below, with the version used by the authors described in parentheses. Please install all packages for command-line usage before proceeding.

Cutadapt (v4.8) https://cutadapt.readthedocs.io/en/stable/ 

BBMap (v38.96) (for clumpify) https://github.com/BioInfoTools/BBMap
**download and place BBMap files into DMS-TRAM-seq/BBMap

STAR (v2.7.1a) https://github.com/alexdobin/STAR

samtools/bcftools (v1.11) https://www.htslib.org/

bedtools (v2.29.2) https://bedtools.readthedocs.io/en/latest/

bedops (v2.4.37) https://bedops.readthedocs.io/en/latest/

python (v3.8)

## Necessary reference data

Though several genomic annotation files are included in this repository that are used in the analysis, two large reference datasets are also necessary (the STAR reference genome for hg38 and chromosome-level fasta files for the hg38 genome sequence). You can download these from zenodo.org/XXXXXXXXX, and they were made as described below. Once the zip files have been downloaded from zenodo:

  ```bash
  unzip STAR_hg38.zip
  mv STAR_hg38 reference_annot/
  unzip fasta_hg38.zip
  mv fasta_hg38 reference_annot/
  ```

Alternatively, you may generate these files yourself as described below:

First, a STAR reference genome for hg38 will need to be built in the folder DMS-TRAM-seq/STAR_hg38. Please follow instructions in the STAR manual at https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf and the reference geome fasta and gtf files can be downloaded from https://www.gencodegenes.org/human/. Only the primary chromosomes are analyzed, so any hg38 version will suffice. Should you already have this STAR reference hg38 built, you can change the "refpath" variable on line 24 of TrimMapDedupProcess.sh to direct there, instead. 

Second, reference fasta files for each hg38 chromosome will need to be placed in DMS-TRAM-seq/reference_annot/fasta, named chrN.fa.

## Installation of this software
  
  ```bash
  git clone https://github.com/farenhem/DMS-TRAM-seq
  cd DMS-TRAM-seq
  ```

# Running the pipeline

## Processing and aligning the sequencing reads

First, the FASTQ files must be appropriately named within the FASTQ folder. They should take the format SampleName_1.fastq.gz and SampleName_2.fastq.gz for read 1 and 2 of SampleName, respectively.

Then, for each sample, run:

  ```bash
  sbatch TrimMapDedupProcess.sh SampleName
  ```
This will fully process the reads, including quality and adapter trimming (cutadapt), sequence-level deduplication (clumpify), alignment (STAR), and filtering the resulting BAM files so that only uniquely-aligned reads are retained in the final "primary" BAM file. The files are now ready for processing into mutational profiles.

## Generating transcriptome-wide mutational profiles

This pipeline for processing BAM files for a set of DMS Map-seq experiments into files of per-region or per-transcript mismatch rates and summary statistics is implemented in 8 sequential steps, the last of which can be repeated to interrogate different elements of the transcriptome.

Each step of the pipeline should successfully run to completion before beginning the subsequent step.

The steps are:

- Calculating mismatch rates for all reference bases:
    The most computationally demanding step, where base mismatch rates calculated transcriptome-wide separately for each chromosome and each sample.  Results  are reported in *.pileup files. 
    ```bash
    ./GWSA1.sh
    ```

- A coverage filter is applied (default: 100 reads per base). Bases are retained after this step only if they pass a filter on depth, which can be set within the shell script.  Output is in bedGraph format.
  ```bash
  ./GWSA2.sh
  ```

- Mismatch rates from DMS-treated samples are "normalized" using non-DMS-treated controls.
  ```bash
  ./GWSA3.sh
  ```

- The per-chromosome, per-sample files from the preceding steps are combined into a single per-sample file.
  ```bash
  ./GWSA4.sh
  ```

- This step of the pipeline, together with the following step, ensures that sites in the transcriptome are retained for further analysis only if they exhibit adequate coverage (i.e. above the specified threshold) for ALL samples/conditions.
  ```bash
  ./GWSA5.sh
  ```

- Together with the previous step, ensure consistently high coverage among regions that are to be compared across conditions.
  ```bash
  ./GWSA6.sh
  ```

- An optional step (NOT recommended, can be skipped entirely) to assess the per-sample distribution of normalized mismatch rates.  These distributions can be used to set a upper bound for a per-sample (as opposed to per-transcript or per-region) Winsorization.
  ```bash
  ./GWSA7.sh
  ```

For now only "canonical" transcripts are considered (i.e. a single annotation
per transcript) .All of the steps below calculate and report summary statistics
for different elements of the transcriptome:
./GWSA8_3utr_assemble.sh # Spliced 3' UTRs.
./GWSA8_5utr_assemble.sh # Spliced 5' UTRs.
./GWSA8_cds_assemble.sh # Spliced coding sequences (i.e. excluding UTRs)
./GWSA8_cds_unspliced.sh # Individual exons (unspliced CDS, UTRs excluded)
./GWSA8_transcript_assemble.sh # Complete spliced transcripts, including UTRs.
