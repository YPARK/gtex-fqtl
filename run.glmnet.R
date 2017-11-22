#!/usr/bin/env Rscript

## Run spike-slab QTL model (tissue by tissue)
argv <- commandArgs(trailingOnly = TRUE) 

if(length(argv) < 4){
    q()
}

library(glmnet)
library(Matrix)
library(methods)
library(fqtl)
library(dplyr)
library(reshape2)

mat.melt <- function(mat, v.name = 'lodds') {
    melt(mat, varnames = c('row', 'col'), value.name = v.name)
}

source('Util.R')
source('Sim.R')

plink.hdr    <- argv[1]  # e.g., 'scratch/temp'
y.file       <- argv[2]  # e.g., 'scratch/temp.y.txt.gz'
snp.out.file <-  argv[3] # e.g., 
pve.out.file <-  argv[4] # e.g., 

run.glmnet <- function(y, x, alpha = 1){
    valid <- !is.na(y)
    xx <- x[valid,,drop=FALSE]
    yy <- as.matrix(y[valid])
    cv.out <- cv.glmnet(x=xx, y=yy, alpha=alpha, nfolds=5)
    ret <- glmnet(x=xx, y=yy, alpha=alpha, lambda=cv.out$lambda.min)$beta
    cat(mean(abs(ret) > 0), '\n')
    return(ret)
}

plink <- read.plink(plink.hdr)

X <- scale(plink$BED, center = TRUE, scale = FALSE)
Y <- scale(as.matrix(read.table(y.file)))

snp.out <- NULL
tis.out <- NULL
n.tis <- dim(Y)[2]

for(aa in c(1, 0.75, 0.5, 0.25)) {
    .beta <- do.call(cbind, apply(Y, 2, run.glmnet, x = X, alpha = aa))
    .pve <- sapply(1:n.tis, function(j) get.pve(Y[, j], .beta[, j], X))

    .tis.out <- data.frame(tis = 1:n.tis, pve = .pve, method = 'glmnet', alpha = aa)

    ## report best tissue effect sizes
    colnames(.beta) <- 1:n.tis
    rownames(.beta) <- plink$BIM[, 2]

    .snp.out <- mat.melt(as.matrix(.beta), 'beta') %>%
        group_by(row) %>% slice(which.max(abs(beta))) %>%
            mutate(alpha = aa)

    snp.out <- rbind(snp.out, as.data.frame(.snp.out))
    tis.out <- rbind(tis.out, .tis.out)

    log.msg('Done alpha = %.2f\n', aa)
}

write.tsv(tis.out, file = gzfile(pve.out.file))
write.tsv(snp.out, file = gzfile(snp.out.file))
