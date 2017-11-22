#!/usr/bin/env Rscript

## Run spike-slab QTL model (tissue by tissue)
argv <- commandArgs(trailingOnly = TRUE) 

if(length(argv) < 4){
    q()
}

plink.hdr    <- argv[1]  # e.g., plink.hdr = 'scratch/simulation/data/13125/3/0.2/3/10/1001-temp'
me.file       <- argv[2]  # e.g., me.file = 'temp.me.txt'
snp.out.file <-  argv[3] # e.g., 
pve.out.file <-  argv[4] # e.g., 


options(stringsAsFactors = FALSE)

library(zqtl)
library(dplyr)
library(reshape2)
source('Sim.R')
source('Util.R')

plink <- read.plink(plink.hdr)
n.ind <- dim(plink$BED)[1]

snps <- plink$BIM[, 2]

me.tab <- read.table(me.file, header = TRUE)
tissues <- sort(unique(me.tab[, 'gene']))
n.snps <- length(snps)
n.tis <- length(unique(tissues))

log.msg('Read %d x %d matrix QTL results\n\n', n.snps, n.tis)

mat.idx <- cbind(match(me.tab[, 'SNP'], snps),
                 match(me.tab[, 'gene'], tissues))

t.mat <- matrix(NA, nrow = n.snps, ncol = n.tis)
beta.mat <- matrix(NA, nrow = n.snps, ncol = n.tis)
se.mat <- matrix(NA, nrow = n.snps, ncol = n.tis)

t.mat[mat.idx] <- me.tab[, 't.stat']
beta.mat[mat.idx] <- me.tab[, 'beta']
se.mat[mat.idx] <- me.tab[, 'beta'] / me.tab[, 't.stat']

mat.melt <- function(mat, v.name = 'lodds') {
    melt(mat, varnames = c('row', 'col'), value.name = v.name)
}

tis.out <- NULL
snp.out <- NULL
rank.out <- NULL

for(k in 1:5) {

    ## slower learning rate & do standardize
    vb.opt <- list(vbiter = 2500, pi.ub = -0, pi.lb = -4, gammax = 1e4,
                   svd.init = TRUE, jitter = 0.1, eigen.tol = 0.01,
                   rate = 1e-2, decay = -0.1, tol = 1e-8, k = k,
                   random.effect = FALSE, do.stdize = TRUE)

    zqtl.out <- fit.zqtl(effect = beta.mat, effect.se = se.mat,
                         n = 0, X = plink$BED, factored = TRUE,
                         options = vb.opt)

    theta <- zqtl.out$param.left$theta %*% t(zqtl.out$param.right$theta)
    theta <- theta * se.mat

    pve <- sapply(1:n.tis, function(j) get.pve.no.y(theta[, j], scale(plink$BED)))

    temp <- data.frame(tis = 1:n.tis,
                       pve = signif(pve, 2),
                       lodds = apply(zqtl.out$param.right$lodds, 1, max),
                       method = 'zfqtl', k)
    
    tis.out <- rbind(tis.out, temp)

    temp <- mat.melt(zqtl.out$param.left$lodds, 'lodds') %>%
        left_join(mat.melt(zqtl.out$param.left$theta, 'theta'), by = c('row', 'col')) %>%
            left_join(mat.melt(zqtl.out$param.left$theta.var, 'theta.var'), by = c('row', 'col')) %>%
                group_by(row) %>% slice(which.max(lodds)) %>%
                    mutate(theta.se = sqrt(theta.var)) %>%
                        select(-theta.var)
    
    snp.out <- rbind(snp.out, data.frame(temp, k = k))

    temp <- data.frame(max.lodds = apply(zqtl.out$param.right$lodds, 2, max),
                       method = 'zfqtl', k)

    rank.out <- rbind(rank.out, temp)
}

write.tsv(tis.out, file = gzfile(pve.out.file))
write.tsv(snp.out, file = gzfile(snp.out.file))

.file.hdr <- glue(dirname(pve.out.file), '/', gsub('.pve.gz', '', basename(pve.out.file)))
rank.out.file <- glue(.file.hdr, '.rank.txt.gz')

write.tsv(rank.out, file = gzfile(rank.out.file))
