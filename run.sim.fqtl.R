#!/usr/bin/env Rscript

## Run spike-slab QTL model (tissue by tissue)
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
library(dplyr)
library(reshape2)
source('Sim.R')
source('Util.R')

mat.melt <- function(mat, v.name = 'lodds') {
    melt(mat, varnames = c('row', 'col'), value.name = v.name)
}

plink <- read.plink(plink.hdr)

X <- scale(plink$BED, center = TRUE, scale = FALSE)
Y <- scale(as.matrix(read.table(y.file)))
n.tis <- dim(Y)[2]

tis.out <- NULL
snp.out <- NULL

for(k in 1:5) {

    vb.opt <- list(vbiter = 2500, pi.ub = -0, pi.lb = -4, gammax = 1e4,
                   svd.init = TRUE, jitter = 0.1,
                   rate = 0.01, decay = -0.1, tol = 1e-6,
                   out.residual = FALSE, k = k)

    fqtl.out <- fqtl.regress(y = Y, x.mean = X, factored = TRUE, options = vb.opt)

    theta <- fqtl.out$mean.left$theta %*% t(fqtl.out$mean.right$theta)

    pve <- sapply(1:n.tis, function(j) get.pve(Y[, j], theta[, j], X))

    temp <- data.frame(tis = 1:n.tis,
                       pve = signif(pve, 2),
                       lodds = apply(fqtl.out$mean.right$lodds, 1, max),
                       method = 'fqtl', k)
    
    tis.out <- rbind(tis.out, temp)

    temp <- mat.melt(fqtl.out$mean.left$lodds, 'lodds') %>%
        left_join(mat.melt(fqtl.out$mean.left$theta, 'theta'), by = c('row', 'col')) %>%
            left_join(mat.melt(fqtl.out$mean.left$theta.var, 'theta.var'), by = c('row', 'col')) %>%
                group_by(row) %>% slice(which.max(lodds)) %>%
                    mutate(theta.se = sqrt(theta.var)) %>%
                        select(-theta.var)
    
    snp.out <- rbind(snp.out, data.frame(temp, k = k))

}

write.tsv(tis.out, file = gzfile(pve.out.file))
write.tsv(snp.out, file = gzfile(snp.out.file))
