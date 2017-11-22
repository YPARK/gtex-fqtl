library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(directlabels)
require(scales)
library(cba)
library(reshape2)

source('figure.util.R')
source('Util.R')
options(stringsAsFactors = FALSE)

k.max <- 50
pip.cutoff <- .95

merge.tabs <- function(...) {
    do.call(rbind, lapply(..., read.table))
}

str.split <- function(...) {
    strsplit(..., split = '[|]')
}

rand.sign <- function() {
    sample(c(-1,1), 1)
}

data.tab <- merge.tabs(glue('result/stat/chr', 1:22, '/', k.max, '/combined.txt.gz'))

colnames(data.tab) <- c('ensg', 'chr', 'tss',
                        'tis.idx', 'tis.names',
                        'tis.theta', 'tis.se', 'tis.lodds',
                        'snps', 'snp.theta', 'snp.se', 'snp.lodds',
                        'k', 'pip')

data.tab <- data.tab %>% filter(pip == pip.cutoff)

## Read tissue colors
tis.colors <- read.table('data/tissues.colors.txt', comment.char = '%')
colnames(tis.colors) <- c('tis.names', 'tis.colors')
tis.colors <- tis.colors %>% filter(tis.names %in% unique(data.tab$tis.names))

## Read expressed gene counts
tis.names <- tis.colors[, 1]
tis.gene.files <- glue('data/expr/',
                       tis.names, '_Analysis.genes.txt.gz')

read.genes <- function(tt) data.frame(read.table(tis.gene.files[tt], header = TRUE), tis.names = tis.names[tt])

genes.tab <- do.call(rbind, lapply(1:length(tis.names), read.genes))
colnames(genes.tab) <- c('ensg', 'tis.names')


################################################################
## Fig 3(a) y = number of genes; x = number of factors (PIP > .95)
tot.count.tab <- genes.tab %>% group_by(tis.names) %>%
    summarize(num.genes.tot = length(unique(ensg))) %>%
    as.data.frame()

n.ind.tab <- read.table('data/individuals/tissues.ind.txt')
colnames(n.ind.tab) <- c('ind', 'pos', 'tis.names')
n.ind.tab <- n.ind.tab %>%
    mutate(tis.names = sapply(tis.names, gsub, pattern = '_Analysis', replacement = '')) %>%
    group_by(tis.names) %>%
    summarize(num.ind = length(unique(ind)))

gene.tis.tab <- data.tab %>%
    select(ensg, tis.idx, tis.names, k) %>%
    mutate(tis.idx = str.split(tis.idx), tis.names = str.split(tis.names)) %>%
    unnest()

gene.count.tab <- gene.tis.tab %>% group_by(tis.names) %>%
    summarize(num.genes = length(unique(ensg))) %>%
    left_join(tot.count.tab, by = 'tis.names')

plt.ngenes.df <- gene.count.tab %>%
    left_join(n.ind.tab, by = 'tis.names') %>%
    mutate(pr.genes = num.genes / num.genes.tot)

plt.ngenes.df$tis.names <- factor(plt.ngenes.df$tis.names, tis.colors[, 1])

dl.method <- list(dl.trans(x = x + .1), 'last.points',
               cex = .5, hjust = 0, vjust = 0, rot = 0)

## fit lm

plt.ngenes.df <-  plt.ngenes.df %>%
    mutate(is.brain =
           sapply(tis.names, substring, first=1, last=5) == 'Brain')

p1 <-
    ggplot(plt.ngenes.df) + theme_bw() +
    geom_smooth(data = plt.ngenes.df,
                aes(x = num.ind, y = num.genes),
                method = 'lm', level = .99, size = .25,
                fill = 'gray80', color = 'gray40', alpha = .5) +
    geom_smooth(data = plt.ngenes.df %>% filter(!is.brain),
                aes(x = num.ind, y = num.genes),
                method = 'lm', level = .99, size = .25,
                fill = '#ccffcc', color = 'green', alpha = .5) +
    geom_point(aes(x = num.ind, y = num.genes, fill = tis.names), pch = 21, show.legend = FALSE) +
    geom_dl(aes(x = num.ind, y = num.genes, group = tis.names, label = tis.names),
            method = dl.method)

p1 <- p1 +
    scale_y_continuous(sec.axis = sec_axis(~ . / max(plt.ngenes.df$num.genes.tot),
                       name = '# eGenes / # expressed'))

p1 <- p1 + scale_fill_manual(values = tis.colors[, 2]) +
    xlab('Sample size') + ylab('# eGenes')

## Fig 3(b) y = number of genes / tissue; x = sample size
factor.per.gene <- data.tab %>%
    filter(pip == pip.cutoff) %>%
    group_by(ensg) %>%
    summarize(num.factors = length(k))

max.k.val <- max(factor.per.gene$num.factors)

factor.per.gene %>%
    filter(num.factors == max.k.val) %>%
    left_join(genes.tab, by = 'ensg')

f.count.tab <- factor.per.gene %>%
    group_by(num.factors) %>%
    summarize(count = length(num.factors))

brk.y <- sort(unique(c(range(f.count.tab$count),
                     c(1, 100, 500, 1000, 2500, 5000))))

brk.x <- sort(unique(c(range(f.count.tab$num.factors),
                       seq(0, 50, 5))))

p2 <- ggplot(f.count.tab, aes(x = num.factors, y = count)) +
    theme_classic() +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = brk.x) +
    scale_y_continuous(breaks = brk.y) +
    xlab('# factors per gene') +
    ylab('# genes')

## Fig 3(c) y = number of genes (or factors); x = number of SNPs
deg.gene.tab <-
    data.tab %>%
    select(ensg, k, tis.names, snps) %>%
    rename(tissues = tis.names) %>%
    mutate(tissues = str.split(tissues)) %>%
    unnest() %>%
    mutate(snps = str.split(snps)) %>%
    unnest() %>% group_by(ensg) %>%
    summarize(snp = length(unique(snps)),
              tissue = length(unique(tissues)))

deg.gene.tab <- deg.gene.tab %>% mutate(unit = 'gene')

deg.factor.tab <-
    data.tab %>%
    select(ensg, k, tis.names, snps) %>%
    mutate(ensg = paste(ensg, k, sep = '-')) %>%
    rename(tissues = tis.names) %>%
    mutate(tissues = str.split(tissues)) %>%
    unnest() %>%
    mutate(snps = str.split(snps)) %>%
    unnest() %>% group_by(ensg) %>%
    summarize(snp = length(unique(snps)),
              tissue = length(unique(tissues)))

deg.factor.tab <- deg.factor.tab %>% mutate(unit = 'factor')

deg.tab <- rbind(deg.gene.tab, deg.factor.tab)

df.snp <- deg.tab %>% group_by(snp, unit) %>%
    summarize(count = length(unique(ensg)))

deg.plot <- function(.df, x.aes) {

    .rng <- range(.df$count)
    .seq <- round(seq(.rng[1], .rng[2], length.out = 10))

    brk <- sort(unique(c(.rng, .seq)))

    ggplot(.df, x.aes) + theme_classic() +
        geom_line(aes(y = count, color = unit)) +
        geom_point(aes(y = count, shape = unit, color = unit)) +
        scale_y_continuous(breaks = brk)
}

deg.prop <- function(stuff) {
    ret <- stuff + scale_shape_manual(values = c(1, 3)) +
        scale_color_manual(values = c('#339933', 'gray40'))

    ret <- ret +
        ylab('# genes (or factors)')

    ret <- ret +
        theme(legend.position = c(.5,1),
              legend.justification = c(.5, .9),
              legend.title = element_blank(),
              legend.background = element_blank())
}

p3 <- deg.prop(deg.plot(df.snp, aes(x = snp))) + xlab('# SNPs')

################################################################
## Fig 3(d) y = number of genes (or factors); x = number of tissues

df.tissue <- deg.tab %>% group_by(tissue, unit) %>%
    summarize(count = length(unique(ensg)))

p4 <- deg.prop(deg.plot(df.tissue, aes(x = tissue))) +
    xlab('# tissues')


################################################################
## Fig 3(e) correlation between factors
.names <- glue('result/stat/chr', 1:22, '/', k.max, '/s2.txt.gz')
s2.tab <- merge.tabs(.names)

colnames(s2.tab) <- c('k1', 'k2', 'z', 'num', 'denom',
                      's2', 'var', 'tr', 'p', 'pip')

s2.tab <- s2.tab %>% filter(k1 != k2)
    
## Just take 99% quantiles and compare
s2.filter <- s2.tab %>%
    filter(num > quantile(num, 0.005), num < quantile(num, 0.995))

p5 <-
    ggplot(s2.filter, aes(x = num / sqrt(s2))) +
    theme_classic() +
    geom_density(fill = 'gray60', color = 'gray60') +
    xlab('correlation z')


g.top <- ggplotGrob(p1)
g.bottom <- match.heights(list(p2, p3, p4, p5))

lay <- rbind(c(1, 1, 1, 1),
             c(1, 1, 1, 1),
             c(2, 3, 4, 5))

pdf('result/figures/fig-gtex-stat.pdf', width = 8, height = 8, useDingbats = FALSE)
grid.arrange(grobs = c(list(g.top), g.bottom),
             layout_matrix = lay)
dev.off()
