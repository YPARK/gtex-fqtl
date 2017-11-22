#!/bin/bash -l

if [ $# -lt 4 ]; then
cat << EOF

plink=\$1
yfile=\$2
output_snp=\$3
output_tis=\$4
tempdir=\$5

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

plink=scratch/simulation/data/29747/1/0.1/1/10/1001-temp
yfile=scratch/simulation/data/29747/1/0.1/1/10/1001-temp.yfull.txt.gz
output_tis=temp.tis.gz

if [ -z $tempdir ]; then
    outdir=$(dirname $output_tis)
    [ -d $outdir ] || mkdir -p $outdir
    tempdir=$(mktemp -d "$outdir/temp.XXXXXX")
fi

[ -d $tempdir ] || mkdir -p $tempdir


# 1. convert plink to vcf
[ -f $tempdir/bma.geno.list ] && rm $tempdir/bma.geno.list

nind=$(cat $plink.fam | wc -l)
ntis=$(zcat $yfile | head -n1| wc -w)
vcf=$tempdir/bma.geno.vcf

[ -f $vcf ] || ./bin/plink --bfile $plink --recode vcf --out $tempdir/bma.geno

for((tis=1; tis<=${ntis}; ++tis)); do
    printf "tis%d\t%s\n" ${tis} ${vcf} >> $tempdir/bma.geno.list
done

# 2. expr list
[ -f $tempdir/bma.exp.list ] && rm $tempdir/bma.exp.list

ensg=$(cat $plink.gene | awk '{ print $1 }')

for((tis=1; tis<=${ntis}; ++tis)); do
    cat $plink.fam | awk '{ if(NR > 1) printf "\t"; printf $1; } END { printf "\n" }' | \
	gzip > $tempdir/bma.exp.${tis}.expr.gz

    zcat $yfile | \
	awk -vtt=${tis} -vG=${ensg} -F'\t' 'BEGIN { printf G; } { printf FS $(tt) } END { printf "\n" }' | \
	gzip >> $tempdir/bma.exp.${tis}.expr.gz

    printf "tis%d\t%s\n" ${tis} $tempdir/bma.exp.${tis}.expr.gz >> $tempdir/bma.exp.list
done

cat $plink.gene | awk -F'\t' '{ print $2 FS ($3 - 1) FS $3 FS $1 }' > $tempdir/bma.gcoord.bed


R --vanilla <<EOF
source('Util.R')
source('bin/utils_eqtlbma.R')
write.tsv(makeGrid("general"), file = gzfile('$tempdir/bma.gridL'))
write.tsv(makeGrid("configs"), file = gzfile('$tempdir/bma.gridS'))
EOF

./bin/eqtlbma_bf --geno $tempdir/bma.geno.list \
		 --gcoord $tempdir/bma.gcoord.bed \
		 --exp $tempdir/bma.exp.list \
		 --cis 1000000 --analys join \
		 --gridL $tempdir/bma.gridL --gridS $tempdir/bma.gridS \
		 --error mvlr --bfs all --out $tempdir/bma.bmabf -v3

