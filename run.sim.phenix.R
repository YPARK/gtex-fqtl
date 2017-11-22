#!/usr/bin/env Rscript

## Run phenix model then call glmnet
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 4){
    q()
}

plink.hdr    <- argv[1]  # e.g., 'scratch/temp'
y.file       <- argv[2]  # e.g., 'scratch/temp.y.txt.gz'
snp.out.file <-  argv[3] # e.g.,
pve.out.file <-  argv[4] # e.g.,

options(stringsAsFactors = FALSE)

library(fqtl)
library(phenix)
source('Sim.R')
source('Util.R')

library(reshape2)

mat.melt <- function(mat, v.name = 'lodds') {
    melt(mat, varnames = c('row', 'col'), value.name = v.name)
}

plink <- read.plink(plink.hdr)

X <- scale(plink$BED, center = TRUE, scale = FALSE)
Y <- scale(as.matrix(read.table(y.file)))
n.tis <- dim(Y)[2]
n.snp <- dim(X)[2]
n.ind <- dim(X)[1]

svd.out <- svd(X)
.idx <- which(svd.out$d > 1e-4)
U.d <- sweep(svd.out$u[, .idx, drop = FALSE], 2, svd.out$d[.idx], `*`)

K <- U.d %*% t(U.d) / (n.snp - 1)

phenix.out <- phenix(Y, K)

Y <- phenix.out$imp

tis.out <- data.frame(tis = 1:n.tis, h2 = phenix.out$h2)

library(glmnet)
library(Matrix)
library(methods)
library(dplyr)

run.glmnet <- function(y, x, alpha = 1){
    valid <- !is.na(y)
    xx <- x[valid,,drop=FALSE]
    yy <- as.matrix(y[valid])
    cv.out <- cv.glmnet(x=xx, y=yy, alpha=alpha, nfolds=5)
    ret <- glmnet(x=xx, y=yy, alpha=alpha, lambda=cv.out$lambda.min)$beta
    cat(mean(abs(ret) > 0), '\n')
    return(ret)
}

snp.out <- NULL

for(aa in c(1, 0.75, 0.5, 0.25)) {

    .beta <- do.call(cbind, apply(Y, 2, run.glmnet, x = X, alpha = aa))
    colnames(.beta) <- 1:n.tis
    rownames(.beta) <- plink$BIM[, 2]

    ## report best tissue effect sizes
    .snp.out <- mat.melt(as.matrix(.beta), 'beta') %>%
        group_by(row) %>% slice(which.max(abs(beta))) %>%
            mutate(alpha = aa)

    snp.out <- rbind(snp.out, as.data.frame(.snp.out))

    log.msg('Done alpha = %.2f\n', aa)
}

write.tsv(tis.out, file = gzfile(pve.out.file))
write.tsv(snp.out, file = gzfile(snp.out.file))
