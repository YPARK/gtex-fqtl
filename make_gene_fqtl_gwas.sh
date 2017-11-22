#!/bin/bash -l

source /broad/software/scripts/useuse > /dev/null
reuse -q BEDTools > /dev/null

GWAS_BED=$1 # e.g., IGAP/chr1.txt.gz
QTL_BIM=$2  # e.g., scratch/data/334/plink.bim

cat $QTL_BIM | awk -F'\t' '{ print "chr" $1 FS ($4 - 1) FS $4 FS $2 FS $5 FS $6 }' | \
    bedtools intersect -a $GWAS_BED -b stdin -wa
