#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 2) q()

k.max <- as.integer(argv[1]) # e.g., 10
out.file <- argv[2]          # e.g., 'result/fig_fdr.pdf'

library(qvalue)
library(ggplot2)
library(directlabels)
library(dplyr)
library(gridExtra)
source('Util.R')
source('figure.util.R')
options(stringsAsFactors = FALSE)
merge.tabs <- function(...) {
    do.call(rbind, lapply(..., read.table))
}

## Take PIP statistics under the null and alternative and estimate FDR
## using Storey's q-value

null.tab <- merge.tabs(glue('result/perm/chr',1:22,'/',k.max,'/snp-null.txt.gz'))
alt.tab <- merge.tabs(glue('result/stat/chr',1:22,'/',k.max,'/snp-null.txt.gz'))

colnames(null.tab) <- c('k', 'lodds')
colnames(alt.tab) <- c('k', 'lodds')

pv <- empPvals(stat = alt.tab[, 2], stat0 = null.tab[, 2])
qv <- qvalue(pv)

## Fig: Show two histograms of PIP values
sigmoid <- function(x) 1/(1+exp(-x))

.df1 <- rbind(null.tab %>% mutate(pip.100 = round(100 * sigmoid(lodds)), stat = 'permuted'),
             alt.tab %>% mutate(pip.100 = round(100 * sigmoid(lodds)), stat = 'observed')) %>%
    group_by(pip.100, stat) %>%
    summarize(count = length(pip.100))

p1 <- ggplot(.df1, aes(x = pip.100, y = count, fill = stat)) +
    theme_classic() +
    geom_vline(xintercept = 50, color = 'orange', lty = 'dashed') +
    geom_vline(xintercept = 95, color = 'orange', lty = 'dashed') +
    geom_bar(stat = 'identity') + xlab('max PIP (%)') +
    coord_cartesian(ylim = c(0, 50000)) +
    theme(legend.title = element_blank(),
          legend.position = c(0.5, 0.6)) +
    scale_x_discrete(limits = seq(0, 100, 10), expand = c(0.01, 0.01))

## Fig: Show q-values for each cutoff

.df2 <- rbind(data.frame(alt.tab, error = pv, stat = 'p-value') %>%
              mutate(pip.100 = round(100 * sigmoid(lodds))),
              data.frame(alt.tab, error = qv$qvalue, stat = 'q-value') %>%
              mutate(pip.100 = round(100 * sigmoid(lodds)))) %>%
    group_by(pip.100, stat) %>%
    summarize(count = length(pip.100), error = max(error)) %>%
    mutate(log.error = -log10(error))

p2 <-
    ggplot(.df2, aes(x = pip.100, y = error, group = stat,
                     shape = stat, color = stat, lty = stat)) +
    geom_vline(xintercept = 50, color = 'orange', lty = 'dashed') +
    geom_hline(data = .df2 %>% filter(pip.100 == 50), aes(yintercept = error),
               color = 'orange', lty = 'dashed') +
    geom_vline(xintercept = 95, color = 'orange', lty = 'dashed') +
    geom_hline(data = .df2 %>% filter(pip.100 == 95), aes(yintercept = error),
               color = 'orange', lty = 'dashed') +
    theme_classic() +
    geom_point() +
    geom_dl(aes(label = stat), method = 'smart.grid') +
    scale_x_discrete(limits = seq(0, 100, 10), expand = c(0.01, 0.01)) +
    theme(axis.title.x=element_blank(), legend.position = 'none')

p2 <- p2 +
    geom_dl(data = .df2 %>% filter(pip.100 == 50),
            aes(label = signif(error, 2)),
            method = list(dl.trans(y = y + 1),'first.points'))+
    geom_dl(data = .df2 %>% filter(pip.100 == 95),
            aes(label = signif(error, 2)),
            method = list(dl.trans(y = y + 1), 'first.points'))+
    scale_color_manual(values = c('black', 'gray40'),
                       limits = c('q-value', 'p-value'))
    
p0 <- ggplot(data.frame(pv), aes(x=pv)) +
    geom_histogram()

pdf(file = out.file, width = 4, height = 4, useDingbats = FALSE)
grid.vcat(list(p2 + scale_y_log10(), p1, p0), heights = c(2,1,1))
dev.off()
