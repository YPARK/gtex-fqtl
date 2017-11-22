#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)

options(stringsAsFactors = FALSE)

rank.tab <- read.table('result/simulation/rank.txt.gz')

colnames(rank.tab) <- c('model', 'fitted', 'gene', 'true', 'h2', 'snp', 'tis', 'rseed')

rank.tab$model <- factor(rank.tab$model, 1:5, paste('K =',1:5))

tis.n <- sort(unique(rank.tab$tis))
rank.tab$tis <- factor(rank.tab$tis, tis.n, paste('#tis of action =', tis.n))

plt.aes <- aes(x = h2, y = mean.fitted, ymin = mean.fitted - 2 * se.fitted,
               ymax = mean.fitted + 2 * se.fitted, color = model)

plt.func <- function(.nsnp) {
    
    .tab <- rank.tab %>% filter(snp == .nsnp)

    .tab <- .tab %>% na.omit() %>%
        group_by(model, h2, true, tis) %>%
            summarize(mean.fitted = mean(fitted),
                      se.fitted = sd(fitted) / sqrt(length(fitted)))
    
    plt <- ggplot(.tab, plt.aes) +
        geom_hline(aes(yintercept = true), lty = 'dashed') +
        geom_linerange() + geom_line() + theme_bw()

    plt <- plt + xlab('Heritability') + ylab('Estimated # factors') +
    facet_grid(true ~ tis, scales = 'free')

    plt + theme(legend.title = element_blank())
}

n.snps <- unique(rank.tab$snp)

plt.list <- lapply(n.snps, plt.func)

r.dir <- 'result/simulation/figures'

plt.files <- paste(r.dir, '/fig2a-rank-snp', n.snps, '.pdf', sep = '')

sapply(1:length(plt.files),
       function(j, ...) ggsave(plt.files[j], plot = plt.list[[j]], ...),
       width = 8, height = 5, useDingbats = FALSE)

