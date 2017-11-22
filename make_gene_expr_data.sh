#!/bin/bash -l

if [ $# -lt 2 ]; then
    cat <<EOF
geneid=\$1
outdir=\$2
EOF
    exit 1
fi

geneid=$1
outdir=$2

cisk=1000
genefile=data/fqtl.genes.txt
matched=data/individuals/matched.txt

geneinfo=$(zless $genefile | head -n ${geneid} | tail -n1)

printf "$geneinfo\n" > $outdir/data.gene

ensg=$(echo $geneinfo | awk '{ print $2 }')
chr=$(echo $geneinfo | awk '{ print $3 }')
tss=$(echo $geneinfo | awk '{ print $4 }')

# skip non-autosomal chromosomes
re='^[0-9]+$'
if ! [[ $chr =~ $re ]] ; then
    printf "%s is not autosomal chromosome\n" $chr
    exit 1
fi

lb=$(echo $tss | awk -vcis=${cisk}000 '{ print ($1 > cis) ? ($1 - cis) : 0; }')
ub=$(echo $tss | awk -vcis=${cisk}000 '{ print ($1 + cis); }')

# 1. extract gene expression values
exprout=$outdir/expr.txt.gz
valid_tissues=$outdir/tissues.txt.gz
[ -f $exprout ] && rm $exprout
[ -f $valid_tissues ] && rm $valid_tissues

tissues=$(cat data/tissues.txt)
tt=0

for tis in $tissues; do

    tt=$((tt + 1))
    genepos=$(zcat data/expr/${tis}.genes.txt.gz | tail -n+2 | awk -vE=${ensg} '{ if($1 == E) print NR }')

    if ! [ -z $genepos ]; then

	printf "Parsing %s ... " data/expr/${tis}.expr.txt

	exprpos=$(cat $matched | awk -vT=${tis} -F'\t' 'BEGIN { n = 0; } ($3 == T) { if(n++ > 0) printf ","; printf $2 }')

	paste -d'\t' <(cat $matched | awk -vT=${tis} -F'\t' '($3 == T) { print $1 }') \
    	      <(zcat data/expr/${tis}.expr.txt.gz| head -n1 | cut -f 2- | tr '\t' '\n' | \
    		       awk -vROWS=${exprpos} -f util_subset_rows.awk) \
    	      <(zcat data/expr/${tis}.expr.txt.gz | tail -n +2 | \
    		       head -n${genepos} | tail -n1 | cut -f 2- | tr '\t' '\n' | \
    		       awk -vROWS=${exprpos} -f util_subset_rows.awk) | \
    	    awk -vTI=${tt} -vT=${tis} -F'\t' '$1 == $2 { print $1 FS $3 FS T FS TI }' | \
    	    gzip >> $exprout

	printf "%d\t%s\n" ${tt} ${tis} | gzip >> $valid_tissues
	printf "done\n"
    fi

done

[ -f $valid_tissues ] || exit 0

# 2. 450 x 44 gene expression matrix
Y=$outdir/Y.txt.gz
Rscript --vanilla make.gene.expr.data.R $exprout $Y || exit 1

# 3. extract genotype information
./bin/plink --bfile genotypes/chr${chr} --make-bed \
	    --chr $chr --from-bp $lb --to-bp $ub \
	    --out $outdir/plink \
	    || exit 1
