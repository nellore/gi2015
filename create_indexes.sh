#!/usr/bin/env bash
## Creates indexes; use hg19
# path to hg19 fasta; derived from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
REF=/scratch0/langmead-fs1/rgasp/hg19.fa
# where to dump output
OUTPUT=/scratch0/langmead-fs1/rgasp

# where to find STAR; static compile of version 2.4.2a is used here
STAR=/scratch0/langmead-fs1/shared/STARv2.4.2a/STAR
# where to find HISAT2 build executable; version 2.0.0-beta is used
HISAT2BUILD=/scratch0/langmead-fs1/shared/hisat2-2.0.0-beta/hisat2-build

# multithread for STAR
CORES=32

mkdir -p $OUTPUT/star
cd $OUTPUT
$STAR --runMode genomeGenerate --genomeDir $OUTPUT/star --genomeFastaFiles $REF --runThreadN $CORES || { echo 'Problem encountered building STAR index'; exit 1; }
$HISAT2BUILD $REF hisat2genome || { echo 'Problem encountered building HISAT2 index'; exit 1; }