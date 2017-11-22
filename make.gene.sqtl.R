#!/usr/bin/env Rscript
## Fit FQTL

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 4) {
    cat('<usage>\nRscript gene.id do.stdize resid.file out.hdr [permutation seed]\n')
    q()
}

source('Util.R')

gene.id   <- as.integer(argv[1]) # gene.id = 11
do.stdize <- as.logical(argv[2]) # do.stdize = TRUE
resid.file <- argv[3]            # resid.file = 'result/fqtl-std/11/fqtl.resid.txt.gz'
out.hdr   <- argv[4]

do.permutation <- FALSE
rseed <- NULL
if(length(argv) > 4) {
    rseed <- as.integer(argv[6])
    log.msg('Permuted fQTL, random seed = %d\n\n', rseed)
    do.permutation <- TRUE
}

################################################################
dir.create(dirname(out.hdr), recursive = TRUE)

if(file.exists(glue(out.hdr, '.tis.lodds.txt.gz'))){
    log.msg('File already exists')
    q()
}

if(!file.exists(glue('scratch/data/', gene.id, '/plink.bed'))) {
    log.msg('No valid SNP\n')
    q()
}

options(stringsAsFactors = FALSE)
library(fqtl)

tis.tab <- read.table(glue('scratch/data/', gene.id, '/tissues.txt.gz'))
plink <- read.plink(glue('scratch/data/', gene.id, '/plink'))
resid.mat <- as.matrix(read.table(resid.file))

get.pve <- function(y, theta, x) {
    var(x %*% theta, na.rm = TRUE) / var(y, na.rm = TRUE)
}

write.tsv.gz <- function(mat, file.name) {
    write.tsv(signif(mat, 4), gzfile(file.name))
}

write.spike.slab <- function(effect, hdr) {
    write.tsv.gz(effect$lodds, glue(hdr, '.lodds.txt.gz'))
    write.tsv.gz(effect$theta, glue(hdr, '.theta.txt.gz'))
    write.tsv.gz(sqrt(effect$theta.var), glue(hdr, '.se.txt.gz'))
}

individuals <- read.table('data/individuals/geno.ind.txt')

## check the order of plink$FAM
fam.ind <- lapply(plink$FAM[, 1], strsplit, split = '-')
fam.ind <- sapply(fam.ind, function(x) paste(x[[1]][1:2], collapse = '-'))
stopifnot(all(fam.ind == individuals))

vb.opt <- list(gammax = 1e4, pi.ub = 0, pi.lb = -2,
               vbiter = 5000, tol = 1e-8,
               jitter = 0.01,
               rate = 5e-3, decay = -5e-3, verbose = TRUE,
               out.residual = FALSE)

X <- scale(plink$BED, scale = do.stdize)
Y <- scale(resid.mat, scale = do.stdize)

if(do.permutation) {
    set.seed(rseed)
    n <- dim(Y)[1]

    Y <- apply(Y, 2,
               function(y) {
                   ret <- matrix(NA, n, 1)
                   ret.pos <- is.finite(y)
                   y.shuf <- y[ret.pos]
                   y.shuf <- sample(y.shuf)
                   ret[ret.pos, 1] <- y.shuf
                   return(ret)
               })

    log.msg('Permuted breaking tissue-tissue correlation\n')
}

sqtl.out <- fqtl.regress(y = Y, x.mean = X, factored = FALSE, options = vb.opt)

## 3. output results
write.tsv(plink$BIM, file = gzfile(glue(out.hdr, '.snp.info.gz')))
write.tsv(tis.tab, file = gzfile(glue(out.hdr, '.tis.info.gz')))

write.spike.slab(sqtl.out$mean, glue(out.hdr, '.effect'))

