#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

twas.file <- argv[1] # e.g., 'result/stat/stwas/fqtl-std/IGAP/21.stwas.combined.txt.gz'
fqtl.file <- argv[2] # e.g., 'result/stat/chr21/50/combined.txt.gz'
out.file <- argv[3]

library(dplyr)
options(stringsAsFactors = FALSE)

twas.cols <- c('ensg', 'chr', 'tss', 'gwas.name',
               'twas.z', 'twas.effect', 'twas.se', 'twas.p', 'factor')

fqtl.cols <- c('ensg', 'chr', 'tss',
               'tis.idx', 'tis.name', 'tis.effect', 'tis.se', 'tis.lodds',
               'snp.name', 'snp.effect', 'snp.se', 'snp.lodds',
               'factor', 'pip.cutoff')

twas.tab <- read.table(twas.file, sep = '\t', col.names = twas.cols)
fqtl.tab <- read.table(fqtl.file, sep = '\t', col.names = fqtl.cols)

out <- fqtl.tab %>% filter(pip.cutoff >= 0.95) %>%
    left_join(twas.tab, by = c('chr', 'ensg', 'tss', 'factor')) %>%
        na.omit()

write.table(out, file = gzfile(out.file,, compression = 9),
            row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
