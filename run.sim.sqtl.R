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

source('Sim.R')
source('Util.R')

plink <- read.plink(plink.hdr)

X <- scale(plink$BED, center = TRUE, scale = FALSE)
Y <- scale(as.matrix(read.table(y.file)))

vb.opt <- list(vbiter = 2500, pi.ub = -1, pi.lb = -2, gammax = 1e4,
               rate = 0.01, decay = -0.1, tol = 1e-6,
               adam.rate.m = 0.5, adam.rate.v = 0.9,
               out.residual = FALSE, k = 5)

sqtl.out <- fqtl.regress(y = Y, x.mean = X, options = vb.opt)
                         
n.tis <- dim(Y)[2]
Beta <- sqtl.out$mean$theta
pve <- sapply(1:n.tis, function(j) get.pve(Y[, j], Beta[, j], X))

tis.out <- data.frame(1:n.tis, signif(pve, 2), method = 'sqtl')
snp.out <- data.frame(snps = plink$BIM[, 2], signif(sqtl.out$mean$theta, 2), method = 'sqtl')

write.tsv(tis.out, file = gzfile(pve.out.file))
write.tsv(snp.out, file = gzfile(snp.out.file))
