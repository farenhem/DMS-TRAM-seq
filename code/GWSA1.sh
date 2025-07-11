#!/bin/sh

# This script is the top-level control script for genome-wide analysis of
# DMS-MaPseq data in different conditions, starting with BAM files (i.e. mapped
# reads) as input.

# This is step 1 of 8, in which mismatch rates and coverage levels are
# tabulated in a strand-specific way for the entire reference genome, starting
# from BAM files.

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

# Process reads into coverage and mismatch files for transcription from each
# strand.
for experiment in ${ctl1} ${ctl2} ${ctl3} ${DMS1} ${DMS2} ${DMS3} ${sctl1} ${sctl2} ${sctl3} ${sDMS1} ${sDMS2} ${sDMS3}
do
    ./getGenomeRatiosStrand_chr.sh ${experiment}
done

exit
