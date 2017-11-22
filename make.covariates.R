#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

## Rearrange covariate matrices
options(stringsAsFactors = FALSE)
source('Util.R')

matched.tab <- read.table('data/individuals/matched.txt')
colnames(matched.tab) <- c('gtex.ind', 'data.pos', 'tis.name', 'ind.pos')
tissues <- read.table('data/tissues.txt')[, 1]

take.tis.covar <- function(tis) {
    covar.raw <- read.table(glue('data/covariates/', tis, '.covariates.txt.gz'))
    cc <- apply(covar.raw[-1, -1], 2, function(c) as.numeric(c))
    n.c <- dim(cc)[1]
    C <- matrix(NA, nrow = 450, ncol = n.c)
    temp <- subset(matched.tab, tis.name == tis)
    C[temp[, 'ind.pos'], ] <- t(cc[, temp[, 'data.pos']])
    return(scale(C))
}

covar.list <- lapply(tissues, take.tis.covar)

dir.create(argv[1], recursive = TRUE)
files <- glue(argv[1], '/', tissues, '.covariates.txt.gz')

lapply(1:length(tissues), function(ti) write.tsv(covar.list[[ti]], gzfile(files[ti])))

