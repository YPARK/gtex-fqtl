#!/usr/bin/env Rscript 
argv <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = FALSE)
tis.ind <- read.table(argv[1])
geno.ind <- read.table(argv[2])

geno.pos <- match(tis.ind[, 1], geno.ind[, 1])

out <- data.frame(tis.ind, geno.pos)

write.table(out, file = argv[3], quote = FALSE,
            sep = '\t', row.names = FALSE, col.names = FALSE)
