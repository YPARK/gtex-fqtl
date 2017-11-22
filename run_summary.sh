#!/bin/bash -l

yfile=$1
plink=$2
covar=$3
output=$4
tempdir=$5

if [ $# -lt 4 ]; then
cat <<EOF
yfile=\$1
plink=\$2
covar=\$3
output=\$4
tempdir=\$5
EOF
    exit 1
fi

if [ -z $tempdir ]; then
    outdir=$(dirname $output)
    [ -d $outdir ] || mkdir -p $outdir
    tempdir=$(mktemp -d "$outdir/temp.XXXXXX")
else
    [ -d $tempdir ] || mkdir -p $tempdir
fi

# calculate summary statistics using MatrixEQTL
Rscript --vanilla run.matrix.qtl.R $yfile $plink $tempdir/me.txt $covar

# clean up results
cat $tempdir/me.txt | tail -n +2 | \
    awk -F'\t' '{ print $1 FS $2 FS $3 FS ($3/$4) FS $5 }' | \
    gzip > $output

[ -d $tempdir ] && rm -r $tempdir
[ -f ${output} ] || exit 1

echo "successfully copied results: $0 $@ --> $output"
