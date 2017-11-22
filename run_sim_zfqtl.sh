#!/bin/bash -l

if [ $# -lt 4 ]; then
cat << EOF

plink=\$1       # e.g., plink=scratch/simulation/data/10287/1/0.15/3/20/1001-temp
yfile=\$2       # e.g., yfile=scratch/simulation/data/10287/1/0.15/3/20/1001-temp.y.txt.gz
output_snp=\$3
output_tis=\$4
tempdir=\$5     # e.g., tempdir=temp

EOF
    exit 1
fi

plink=$1
yfile=$2
output_snp=$3
output_tis=$4
tempdir=$5

cat << EOF

plink=$1
yfile=$2
output_snp=$3
output_tis=$4
tempdir=$5

EOF

if [ -z $tempdir ]; then
    outdir=$(dirname $output_tis)
    [ -d $outdir ] || mkdir -p $outdir
    tempdir=$(mktemp -d "$outdir/temp.XXXXXX")
fi

[ -d $tempdir ] || mkdir -p $tempdir


# 1. caculate summary statistics

Rscript --vanilla run.matrix.qtl.R $yfile $plink $tempdir/me.txt || exit 1

printf "[%s] Generated: $tempdir/me.txt\n\n\n" "$(date)"

# 2. run z-FQTL

Rscript --vanilla run.sim.zfqtl.R $plink $tempdir/me.txt $output_snp $output_tis || exit 1

printf "[%s] Finished: $0 --> %s and %s\n\n" "$(date)" $output_tis $output_snp
