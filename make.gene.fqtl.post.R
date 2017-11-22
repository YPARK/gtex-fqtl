#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 5) {
    q()
}

fqtl.hdr <- argv[1] # e.g., fqtl.hdr <- 'result/fqtl/334/6/fqtl'
data.dir <- argv[2] # e.g., data.dir <- 'scratch/data/334/'
single.file <- argv[3]
pair.file <- argv[4]
out.file <- argv[5]

if(file.exists(out.file) && file.exists(single.file) && file.exists(pair.file)) {
   q()
}


source('Util.R')
library(fqtl)

plink.hdr <- glue(data.dir, '/plink')
options(stringsAsFactors = FALSE)

gene.info <- read.table(glue(data.dir, '/data.gene'))
gene.info <- gene.info[2:4]

plink <- read.plink(plink.hdr)
X <- scale(plink$BED)
n <- dim(X)[1]
LD <- t(X) %*% X / (n - 1)

snp.lodds <- read.table(glue(fqtl.hdr,'.snp.lodds.txt.gz'))
snp.theta <- as.matrix(read.table(glue(fqtl.hdr,'.snp.theta.txt.gz')))
snp.se <- as.matrix(read.table(glue(fqtl.hdr,'.snp.se.txt.gz')))
snp.info <- as.matrix(read.table(glue(fqtl.hdr,'.snp.info.gz')))

tis.lodds <- read.table(glue(fqtl.hdr,'.tis.lodds.txt.gz'))
tis.theta <- as.matrix(read.table(glue(fqtl.hdr,'.tis.theta.txt.gz')))
tis.se <- as.matrix(read.table(glue(fqtl.hdr,'.tis.se.txt.gz')))
tis.info <- as.matrix(read.table(glue(fqtl.hdr,'.tis.info.gz')))

tis.info[, 1] <- as.integer(tis.info[, 1])
tis.info[, 2] <- sapply(tis.info[, 2], gsub, pattern = '_Analysis', replacement = '')

max.lodds <- apply(snp.lodds, 2, max)

## 1. basic statistics

## stat = theta1' * R * theta2
## stat ~ N(theta1' * R * theta2, Var)
## Var  = theta1' * R * theta1 + theta2' * R * theta2 + p
##     += theta1' * (R * V2) * theta1
##     += theta2' * (R * V1) * theta2
##     += tr[R*V1] + tr[R*V2] + tr[R*V1*R*V2]

take.z.pair.stat <- function(ff, lo.cutoff) {

    .lodds1 <- snp.lodds[, ff[1], drop = FALSE]
    .lodds2 <- snp.lodds[, ff[2], drop = FALSE]

    .snps <- which(.lodds1 > lo.cutoff | .lodds2 > lo.cutoff)
    .p <- length(.snps)

    .theta1 <- snp.theta[.snps, ff[1], drop = FALSE]
    .theta2 <- snp.theta[.snps, ff[2], drop = FALSE]

    .var1 <- diag(as.numeric(snp.se[.snps, ff[1], drop = FALSE]^2), nrow = .p, ncol = .p)
    .var2 <- diag(as.numeric(snp.se[.snps, ff[2], drop = FALSE]^2), nrow = .p, ncol = .p)

    .ld <- LD[.snps, .snps, drop = FALSE]

    .ld1 <- .ld * .var1
    .ld2 <- .ld * .var2
    .ld12 <- .ld1 * .ld2

    .num <- as.numeric(t(.theta1) %*% .ld %*% .theta2)

    .s2 <- t(.theta1) %*% (.ld) %*% (.theta1) +
        t(.theta2) %*% (.ld) %*% (.theta2)
    .s2 <- as.numeric(.s2)

    .var <- t(.theta1) %*% (.ld2 * .ld) %*% (.theta1) +
        t(.theta2) %*% (.ld1 * .ld) %*% (.theta2)
    .var <- as.numeric(.var)

    .tr <- .p + sum(diag(.ld1)) + sum(diag(.ld2)) + sum(diag(.ld12))
    .tr <- as.numeric(.tr)

    .denom <- .s2 + .var + .tr

    .denom <- as.numeric(sqrt(.denom))

    data.frame(f1 = ff[1], f2 = ff[2],
               z = .num / .denom,
               z.num = .num,
               z.denom = .denom,
               s2 = .s2,
               var = .var,
               tr = .tr,
               p = .p)
}




pair.output <- lapply(c(0.1, 0.5, 0.9, 0.95), function(pip.cutoff) {
    lo.cutoff <- log(pip.cutoff) - log(1 - pip.cutoff)
    factors <- which(max.lodds > lo.cutoff)
    ret <- NULL

    if(length(factors) > 1) {
        ret.0 <- do.call(rbind, apply(rbind(factors, factors), 2, take.z.pair.stat, lo.cutoff = lo.cutoff))
        factor.pairs <- combn(factors, 2)
        ret <- do.call(rbind, apply(factor.pairs, 2, take.z.pair.stat, lo.cutoff = lo.cutoff))
        ret <- data.frame(rbind(ret, ret.0), pip.cutoff)
    }

    return(ret)
})

pair.output <- do.call(rbind, pair.output)



single.output <- lapply(c(0.1, 0.5, 0.9, 0.95), function(pip.cutoff) {
    lo.cutoff <- log(pip.cutoff) - log(1 - pip.cutoff)
    data.frame(n.factors = sum(max.lodds > lo.cutoff), pip.cutoff)
})

single.output <- do.call(rbind, single.output)


write.tsv(pair.output, file = gzfile(pair.file))
write.tsv(single.output, file = gzfile(single.file))


## 2. output selected SNPs and tissues
collapse <- function(...) paste(..., collapse = '|')

collapse.result <- function(f, lo.cutoff) {

    .tis <- which(tis.lodds[, f] > lo.cutoff)
    .snp <- which(snp.lodds[, f] > lo.cutoff)

    if(length(.tis) < 1 || length(.snp) < 1) {
        return(NULL)
   } else {
       ret <- data.frame(tis.idx = collapse(tis.info[.tis, 1]),
                         tis.name = collapse(tis.info[.tis, 2]),
                         tis.theta = collapse(signif(tis.theta[.tis, f], 2)),
                         tis.se = collapse(signif(tis.se[.tis, f], 2)),
                         tis.lo = collapse(signif(tis.lodds[.tis, f], 2)),
                         snp.name = collapse(snp.info[.snp, 2]),
                         snp.theta = collapse(signif(snp.theta[.snp, f], 2)),
                         snp.se = collapse(signif(snp.se[.snp, f], 2)),
                         snp.lo = collapse(signif(snp.lodds[.snp, f], 2)),
                         factor = f)
       return(ret)
   }
}

out <- lapply(c(0.5, 0.9, 0.95), function(pip.cutoff) {
    lo.cutoff <- log(pip.cutoff) - log(1 - pip.cutoff)
    factors <- which(max.lodds > lo.cutoff)
    ret <- do.call(rbind, lapply(factors, collapse.result, lo.cutoff = lo.cutoff))
    if(!is.null(ret)){
        ret <- cbind(ret, pip.cutoff)
    }
    return(ret)
})

out <- do.call(rbind, out)
if(!is.null(out)){
    out <- data.frame(ensg = gene.info[1], chr = gene.info[2], tss = gene.info[3], out)
}
write.tsv(out, file = gzfile(out.file))
