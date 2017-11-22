#!/usr/bin/env Rscript
## construct individual by tissue expression matrix
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 2) {
    q()
}

expr.raw.file <- argv[1] # e.g., expr.raw.file = 'scratch/data/10001/expr.txt.gz'
out.file <- argv[2]      # e.g., 'scratch/expr.txt.gz'

library(dplyr)
source('Util.R')

options(stringsAsFactors = FALSE)
matched.tab <- read.table('data/individuals/matched.txt')
colnames(matched.tab) <- c('gtex.ind', 'data.pos', 'tis', 'ind.pos')

expr.raw.dat <- read.table(expr.raw.file)
colnames(expr.raw.dat) <- c('gtex.ind', 'value', 'tis', 'tis.pos')

expr.raw.dat <- left_join(expr.raw.dat, matched.tab, by = 'gtex.ind')

Y <- matrix(NA, 450, 44)
Y[cbind(expr.raw.dat$ind.pos, expr.raw.dat$tis.pos)] <- expr.raw.dat$value

write.tsv(Y, file = gzfile(out.file))

