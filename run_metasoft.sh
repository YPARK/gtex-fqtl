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
file=$2
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

R --vanilla <<EOF
options(stringsAsFactors = FALSE);
me.file = '$tempdir/me.txt'
me.tab = read.table(me.file, header = TRUE)

snps = sort(unique(me.tab[, 'SNP']))
genes = sort(unique(me.tab[, 'gene']))
n.snps = length(unique(snps))
n.genes = length(unique(genes))

beta.idx = cbind(match(me.tab[, 'SNP'], snps),
                 2 * match(me.tab[, 'gene'], genes) - 1)

se.idx = cbind(match(me.tab[, 'SNP'], snps),
                 2 * match(me.tab[, 'gene'], genes))

out = matrix(NA, nrow = n.snps, ncol = n.genes * 2)

out[beta.idx] = me.tab[, 'beta']
out[se.idx] = me.tab[, 'beta'] / me.tab[, 't.stat']

write.table(data.frame(snps, out), file = '$tempdir/input.txt',
row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)

cat(genes, file = '$tempdir/me.traits.txt', sep = '\n')
EOF

# 2. Run metasoft method

java -jar bin/Metasoft.jar \
    -pvalue_table bin/HanEskinPvalueTable.txt \
    -mvalue \
    -mcmc_burnin 10000 \
    -mcmc_sample 50000 \
    -mvalue_method mcmc \
    -mvalue_p_thres 1e-2 \
    -input $tempdir/input.txt \
    -log $tempdir/log.txt \
    -seed 1667 \
    -output $tempdir/output.txt || \
    exit 1

printf "[%s] Successfully run JAVA [%s]\n\n" "$(date)" "$tempdir/output.txt"

# 3. Extract SNP statistics
# Just report RE2 (GTEx default)
cat $tempdir/output.txt | \
    tail -n+2 | \
    awk -F'\t' '{ print $1 FS $9 FS $10 FS $11 }' | \
    gzip > $output_snp

# 4. Choose the best SNP and crop p-values and m-values

best=$(zcat $output_snp | sort -k2g | head -n1)
bestSNP=$(echo $best | awk '{ print $1 }')

cat $tempdir/output.txt | awk -vSNP=$bestSNP -F'\t' '$1 == SNP {
nstudy = $2
rsid = $1
for(j=17; j<(17+nstudy); ++j) {
print rsid FS (j-16) FS $(j) FS $(j + nstudy)
}
}' | gzip > $output_tis

printf "[%s] Finished: $0 --> %s and %s\n\n" "$(date)" $output_tis $output_snp
