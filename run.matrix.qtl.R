#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 3) {
    q()
}

y.file <- argv[1]      # e.g., y.file = 'scratch/simulation/data/13125/3/0.2/3/10/1001-temp.y.txt.gz'
plink.hdr <- argv[2]   # e.g., plink.hdr = 'scratch/simulation/data/13125/3/0.2/3/10/1001-temp'
me.out.file <- argv[3] # e.g., me.out.file = 'temp.me.txt'
covar.file <- NULL

if(length(argv) > 3) {
    covar.file <- argv[4]
}


source('Sim.R')
source('Util.R')
source('Util.matrix.qtl.R')

## run MatrixEQTL
options(stringsAsFactors = FALSE)
library(fqtl)
plink <- read.plink(plink.hdr)
Y <- center(as.matrix(read.table(y.file)))

if(!is.null(covar.file)){
    C <- center(as.matrix(read.table(covar.file)))
} else {
    C <- NULL
}

snp.names <- plink$BIM[, 2]
individuals <- plink$FAM[, 1]

me <- run.matrix.qtl(X = center(plink$BED), Y, C, snp.names, individuals, me.out.file)

## eqtl.tab <- read.table(me$param$output_file_name, header = TRUE)
## system(paste('gzip', me.out.file))
