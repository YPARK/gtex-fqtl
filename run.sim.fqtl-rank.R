#!/usr/bin/env Rscript


## Check fQTL rank
argv <- commandArgs(trailingOnly = TRUE) 

if(length(argv) < 3){
    q()
}

plink.hdr    <- argv[1]  # e.g., 'scratch/simulation/data/13125/1/0.15/3/30/1001-temp'
y.file       <- argv[2]  # e.g., 'scratch/simulation/data/13125/1/0.15/3/30/1001-temp.y.txt.gz'
max.lodds.out.file <-  argv[3] # e.g., 

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

    temp <- data.frame(max.lodds = apply(fqtl.out$mean.right$lodds, 2, max),
                       method = 'fqtl', k)
    
    tis.out <- rbind(tis.out, temp)
}

write.tsv(tis.out, file = gzfile(max.lodds.out.file))
