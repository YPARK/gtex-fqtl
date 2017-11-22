#!/usr/bin/env Rscript
## Calculate column-wise statistics
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 2) {
    q()
}

in.file <- argv[1]
out.file <- argv[2]

library(dplyr)
library(reshape2)
dat <- read.table(in.file, header = TRUE)

out <- melt(dat) %>% group_by(variable) %>%
    summarize(mean = mean(value),
              sd = sd(value),
              q.025 = quantile(value, probs = 0.025),
              q.05 = quantile(value, probs = 0.05),
              q.95 = quantile(value, probs = 0.95),
              q.975 = quantile(value, probs = 0.975))

write.tsv <- function(...) {
    write.table(..., col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
}

write.tsv(as.data.frame(out), file = out.file)
