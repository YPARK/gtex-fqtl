#!/usr/bin/env Rscript
                                        # Generate simulation data
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 6) { q() }

K            <- as.integer(argv[1]) # e.g., K = 1
h2           <- as.numeric(argv[2]) # e.g., 0.1
n.causal.snp <- as.integer(argv[3]) # number of causal SNPs for each rank
n.causal.tis <- as.integer(argv[4]) # number of tissue of action for each rank
rseed        <- as.integer(argv[5]) # e.g., rseed = 1667
plink.hdr    <- argv[6]             # e.g., plink.hdr = 'temp'

library(fqtl)
library(dplyr)
source('Sim.R')
source('Util.R')

options(stringsAsFactors = FALSE)
matched.tab <- read.table('data/individuals/matched.txt')
colnames(matched.tab) <- c('ind', 'data.pos', 'tis', 'ind.pos')
tissues <- unique(matched.tab$tis)
n.tis <- length(tissues)

matched.tab <- matched.tab %>% mutate(tis.pos = match(tis, tissues))

plink <- read.plink(plink.hdr)
snps <- plink$BIM[, 2]
X <- scale(plink$BED, center = TRUE, scale = FALSE)

sim.data <- sample.tissues(X, n.tis, n.causal.snp, n.causal.tis,
                           K, h2, rseed)

## write simulated data
out.hdr <- plink.hdr

y.full.file <- glue(out.hdr, '.yfull.txt.gz')
y.file <- glue(out.hdr, '.y.txt.gz')
eta.file <- glue(out.hdr, '.eta.txt.gz')
causal.file <- glue(out.hdr, '.causal.txt.gz')
settings.file <- glue(out.hdr, '.settings.txt.gz')

Y <- add.missing.values(sim.data$output, matched.tab)
write.tsv(as.data.frame(sim.data$output), file = gzfile(y.full.file))
write.tsv(as.data.frame(Y), file = gzfile(y.file))
write.tsv(as.data.frame(sim.data$eta), file = gzfile(eta.file))
write.tsv(data.frame(sim.data$causal.pairs, snps[sim.data$causal.pairs[, 1]]),
          file = gzfile(causal.file))
write.tsv(data.frame(K, h2, n.causal.snp, n.causal.tis, rseed),
          file = gzfile(settings.file))
