#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE) 

if(length(argv) < 2){
    q()
}

options(stringsAsFactors = FALSE)
source('Util.R')
library(dplyr)
sim.tag <- argv[1]         # e.g., sim.tag = '9913/1/0.2/5/30/1001'
out.file <- argv[2]

max.lodds.file <-  glue('scratch/simulation/result/fqtl/', sim.tag, '-fqtl.lodds.gz')
cutoff <- 0
max.lodds.tab <- read.table(max.lodds.file)

colnames(max.lodds.tab) <- c('lodds', 'method', 'model.rank')
out <- max.lodds.tab %>% group_by(model.rank) %>%
    summarize(fit.rank <- sum(lodds > cutoff))

write.tsv(data.frame(out, tag = gsub('/', '\t', sim.tag)),
          file = gzfile(out.file))

