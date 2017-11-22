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

if [ -z $tempdir ]; then
    outdir=$(dirname $output_tis)
    [ -d $outdir ] || mkdir -p $outdir
    tempdir=$(mktemp -d "$outdir/temp.XXXXXX")
fi

[ -d $tempdir ] || mkdir -p $tempdir


# 1. caculate summary statistics

Rscript --vanilla run.matrix.qtl.R $yfile $plink $tempdir/me.txt

printf "[%s] Generated: $tempdir/me.txt\n\n\n" "$(date)"

# 2. Run MASH

R --vanilla << EOF

me.file <- '$tempdir/me.txt'
snp.file <- '$output_snp'
tis.file <- '$output_tis'

library(ashr)
library(mashr)

source('Util.R')

options(stringsAsFactors = FALSE)
me.tab <- read.table(me.file, header = TRUE)

snps <- sort(unique(me.tab[, 'SNP']))
genes <- sort(unique(me.tab[, 'gene']))
n.snps <- length(unique(snps))
n.genes <- length(unique(genes))

log.msg('Read %d x %d matrix QTL results\n\n', n.snps, n.genes)

mat.idx <- cbind(match(me.tab[, 'SNP'], snps),
                 match(me.tab[, 'gene'], genes))

t.mat <- matrix(NA, nrow = n.snps, ncol = n.genes)
beta.mat <- matrix(NA, nrow = n.snps, ncol = n.genes)
se.mat <- matrix(NA, nrow = n.snps, ncol = n.genes)

t.mat[mat.idx] <- me.tab[, 't.stat']
beta.mat[mat.idx] <- me.tab[, 'beta']
se.mat[mat.idx] <- me.tab[, 'beta'] / me.tab[, 't.stat']

mash.data <- set_mash_data(beta.mat, se.mat)
U.c <- cov_canonical(mash.data)

mash.out <- mash(mash.data, U.c, optmethod = 'cxxMixSquarem')

## a. extract SNP statistics

fdr.out <- get_lfdr(mash.out)
fsr.out <- get_lfsr(mash.out)
min.fdr <- apply(fdr.out, 1, min)
min.fsr <- apply(fsr.out, 1, min)

snp.out <- data.frame(snps, min.fdr, min.fsr)

write.table(snp.out, file = gzfile(snp.file),
            quote = FALSE, row.names = FALSE, col.names = FALSE)



## b. choose the best SNP and crop posterior probabilities of genes

best.snp <- which.min(min.fsr)

tis.out <- data.frame(snps[best.snp],
                      genes,
                      fdr.out[best.snp, ],
                      fsr.out[best.snp, ])

write.table(tis.out, file = gzfile(tis.file),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

EOF

printf "[%s] Finished: $0 --> %s and %s\n\n" "$(date)" $output_tis $output_snp
