#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 4) q()

k.max <- as.integer(argv[1])      # e.g., k.max = 50
pip.cutoff <- as.numeric(argv[2]) # e.g., pip.cutoff = 95
stat1.file <- argv[3]             # e.g., 'result/fig_stat1.pdf'
stat2.file <- argv[4]             # e.g., 'result/fig_stat2.pdf'

if(k.max < 2) q()

library(ggplot2)
library(dplyr)

source('Util.R')
options(stringsAsFactors = FALSE)

merge.tabs <- function(...) {
    do.call(rbind, lapply(..., read.table))
}

data.tab <- merge.tabs(glue('result/stat/chr', 1:22, '/', k.max, '/combined.txt.gz'))

colnames(data.tab) <- c('ensg', 'chr', 'tss',
                        'tis.idx', 'tis.names', 'tis.theta', 'tis.se', 'tis.lodds',
                        'snps', 'snp.theta', 'snp.se', 'snp.lodds', 'k', 'pip')

pip.cutoff <- pip.cutoff / 100

## Fig: histogram of #factors

count.tab <- data.tab %>%
    filter(pip == pip.cutoff) %>%
    group_by(ensg) %>%
    summarize(num.factors = length(k))

s1.tab <- count.tab %>%
    group_by(num.factors) %>%
    summarize(count = length(num.factors))

summary(glm(data = s1.tab, formula = count ~ num.factors +1, family = Gamma(link = 'log')))

brk <- sort(unique(c(range(s1.tab$count), c(1, 10, 50, 100, 500, 1000))))

p1 <- ggplot(s1.tab, aes(x = num.factors, y = count)) +
    theme_classic() +
    geom_point() +
    scale_x_continuous(breaks = seq(1,50)) +
    scale_y_log10(breaks = brk) + xlab('number of factors in one gene') +
    theme(axis.text = element_text(size=6))

pdf(file = stat1.file, width = 2.5, height = 3, useDingbats = FALSE)
print(p1)
dev.off()

## Fig: histogram of factor-factor covariance

s2.tab <- merge.tabs(glue('result/stat/chr', 1:22, '/', k.max, '/s2.txt.gz'))

colnames(s2.tab) <- c('k1', 'k2', 'z', 'num', 'denom',
                      's2', 'var', 'tr', 'p', 'pip')

s2.tab <- s2.tab %>% mutate(self = ifelse(k1 == k2, 'within', 'between'))

## Just take 95% quantiles and compare

s2.filter <- s2.tab %>% group_by(self) %>%
    filter(num > quantile(num, 0.025), num < quantile(num, 0.975))

p2 <- ggplot(s2.filter, aes(y = num / sqrt(s2), x = self, color = self, fill = self)) +
    theme_classic() +
    geom_hline(yintercept = 0, color = 'gray60', lty = 'dashed') +
    geom_violin() +
    ylab('z-score') +
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.title.x = element_blank())

p2 <- p2 + scale_color_manual(values = c('orange', 'gray40')) +
    scale_fill_manual(values = c('orange', 'gray40'))

pdf(file = stat2.file, width = 2, height = 3, useDingbats = FALSE)
print(p2)
dev.off()
