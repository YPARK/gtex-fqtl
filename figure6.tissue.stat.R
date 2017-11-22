#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 5) q()

k.max <- as.integer(argv[1])      # e.g., k.max = 50
pip.cutoff <- as.numeric(argv[2]) # e.g., pip.cutoff = 95

num.gene.file <- argv[3]
num.tis.file <- argv[4]
tis.pair.file <- argv[5]


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

merge.tabs <- function(...) {
    do.call(rbind, lapply(..., read.table))
}

str.split <- function(...) {
    strsplit(..., split = '[|]')
}

rand.sign <- function() {
    sample(c(-1,1), 1)
}

pip.cutoff <- pip.cutoff / 100

################################################################
data.tab <- merge.tabs(glue('result/stat/chr', 1:22, '/', k.max, '/combined.txt.gz'))

colnames(data.tab) <- c('ensg', 'chr', 'tss',
                        'tis.idx', 'tis.names', 'tis.theta', 'tis.se', 'tis.lodds',
                        'snps', 'snp.theta', 'snp.se', 'snp.lodds', 'k', 'pip')

data.tab <- data.tab %>% filter(pip == pip.cutoff)

## Read tissue colors
tis.colors <- read.table('data/tissues.colors.txt', comment.char = '%')
colnames(tis.colors) <- c('tis.names', 'tis.colors')
tis.colors <- tis.colors %>% filter(tis.names %in% unique(data.tab$tis.names))

## Read expressed gene counts
tis.names <- tis.colors[, 1]
tis.gene.files <- glue('data/expr/', tis.names, '_Analysis.genes.txt.gz')

read.genes <- function(tt) data.frame(read.table(tis.gene.files[tt], header = TRUE),
                                      tis.names = tis.names[tt])

genes.tab <- do.call(rbind, lapply(1:length(tis.names), read.genes))
colnames(genes.tab) <- c('ensg', 'tis.names')

tot.count.tab <- genes.tab %>% group_by(tis.names) %>%
    summarize(num.genes.tot = length(unique(ensg))) %>%
    as.data.frame()

################################################################
## Fig: histogram of #genes per tissue
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

p1 <- ggplot(plt.ngenes.df) + theme_bw() +
    geom_point(aes(x = num.ind, y = num.genes, fill = tis.names), pch = 21, show.legend = FALSE) +
    geom_dl(aes(x = num.ind, y = num.genes, group = tis.names, label = tis.names),
            method = dl.method)

p1 <- p1 + scale_y_continuous(sec.axis = sec_axis(~ . / max(plt.ngenes.df$num.genes.tot),
                                  name = 'heritable / max expressed genes'))

p1 <- p1 + scale_fill_manual(values = tis.colors[, 2]) +
    xlab('Sample size') + ylab('# of heritable genes')

pdf(file = num.gene.file, width = 8, height = 3, useDingbats = FALSE)
print(p1)
dev.off()



## Fig: #tissues and #snps per factor (or gene)
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


df.p2.tissue <- deg.tab %>% group_by(tissue, unit) %>%
    summarize(count = length(unique(ensg)))

summary(glm(data = df.p2.tissue %>% filter(unit == 'factor') %>% as.data.frame(),
            formula = count ~ tissue + 1, family = Gamma(link = 'log')))

summary(glm(data = df.p2.tissue %>% filter(unit == 'gene') %>% as.data.frame(),
            formula = count ~ tissue + 1, family = Gamma(link = 'log')))


df.p2.snp <- deg.tab %>% group_by(snp, unit) %>%
    summarize(count = length(unique(ensg)))

summary(glm(data = df.p2.snp %>% filter(unit == 'factor') %>% as.data.frame(),
            formula = count ~ snp + 1, family = Gamma(link = 'log')))

summary(glm(data = df.p2.snp %>% filter(unit == 'gene') %>% as.data.frame(),
            formula = count ~ snp + 1, family = Gamma(link = 'log')))



df.p2.xy <- deg.tab %>% group_by(tissue, unit) %>%
    summarize(snp.sd = sd(snp), snp.n = length(snp), snp = mean(snp))

summary(lm(data = df.p2.xy %>% filter(unit == 'factor') %>% as.data.frame(),
            formula = snp ~ tissue + 1))

summary(lm(data = df.p2.xy %>% filter(unit == 'gene') %>% as.data.frame(),
            formula = snp ~ tissue + 1))



p2.tis.stat <-
    deg.tab %>% group_by(unit) %>%
    summarize(med = median(tissue),
              mean = mean(tissue),
              sd = sd(tissue),
              q.25 = quantile(tissue, .25),
              q.75 = quantile(tissue, .75))

p2.snp.stat <-
    deg.tab %>% group_by(unit) %>%
    summarize(med = median(snp),
              mean = mean(snp),
              sd = sd(snp),
              q.25 = quantile(snp, .25),
              q.75 = quantile(snp, .75))

deg.plot <- function(.df, x.aes) {

    brk <- sort(unique(c(range(.df$count), c(1, 10, 50, 100, 500, 1000))))

    ggplot(.df, x.aes) + theme_classic() +
        geom_point(aes(y = count, shape = unit, color = unit)) +
        scale_y_log10(breaks = brk)
}

stat.plot <- function(.stat) {
    aes.h <- aes(y = unit, x = mean, xmin = pmax(mean - sd, 0), xmax = mean + sd, color = unit)
    aes.s <- aes(y = unit, yend = unit, x = q.25, xend = q.75, color = unit)
    .r.data <- .stat %>% select(-sd) %>% melt(id.vars = c('unit'))
    aes.r <- aes(y = unit, x = value, label = signif(value,2), color = unit)

    ret <- ggplot(.stat) + theme_void() + theme(legend.position = 'none') +
        geom_errorbarh(aes.h, height = .1) +
            geom_segment(aes.s, size = 1.5) +
                geom_text_repel(data = .r.data, aes.r, nudge_y = .5, size = 3) +
                    geom_point(aes(y = unit, x = mean, shape = unit, color = unit), size = 2) +
                        geom_dl(aes(x = mean + sd, y = unit, label = unit, color = unit),
                                method = list(dl.trans(x = x + 1), 'last.points', vjust = 0, hjust = 1,
                                    cex = .5))
    return(ret)
}


p2.1 <- deg.plot(df.p2.snp, aes(x = snp))

p2.1.1 <- stat.plot(p2.snp.stat) + scale_x_continuous(limits = c(0, max(df.p2.snp$snp)))

p2.2 <- deg.plot(df.p2.tissue, aes(x = tissue))

p2.2.1 <- stat.plot(p2.tis.stat) + scale_x_continuous(limits = c(0, max(df.p2.tissue$tissue)))

p2.xy.aes <- aes(x = tissue, y = snp,
                 ymin = snp - 2 *snp.sd/sqrt(snp.n),
                 ymax = snp + 2 *snp.sd/sqrt(snp.n),
                 group = unit, shape = unit, color = unit)

p2.xy <- ggplot(df.p2.xy, p2.xy.aes) +
    theme_classic() +
    geom_abline(slope = 1, intercept = 0, color = 'gray80', linetype = 'dotted') +
    geom_linerange() +
    geom_point()

p2.prop <- function(stuff) {
    stuff + scale_shape_manual(values = c(1, 3)) +
        scale_color_manual(values = c('#339933', 'gray40')) +
    theme(legend.position = c(.5,1), legend.justification = c(.5, .9),
          legend.title = element_blank(), legend.background = element_blank())
}


g2.top <- match.heights(lapply(list(p2.1, p2.2, p2.xy), p2.prop))
p2.bottom <- lapply(list(p2.1.1, p2.2.1, ggplot() + theme_void()), p2.prop)
p2.bottom <- lapply(p2.bottom, function(p) p + theme(legend.position = 'none'))
g2.bottom <- match.heights(p2.bottom)

g2.bottom[[1]]$widths[2:5] <- g2.top[[1]]$widths[2:5]
g2.bottom[[2]]$widths[2:5] <- g2.top[[2]]$widths[2:5]

pdf(file = num.tis.file, width = 8, height = 4, useDingbats = FALSE)
grid.arrange(grobs = c(g2.top, g2.bottom), nrow = 2, heights = c(3,1))
dev.off()


################################################################
## Fig: tissue-tissue Jaccard
tis.names <- tis.colors[, 1]
tis.pairs <- combn(tis.names, 2)

take.jaccard <- function(pp) {

    g1 <- gene.tis.tab %>% filter(tis.names == pp[1]) %>%
        select(ensg) %>% unique()

    g2 <- gene.tis.tab %>% filter(tis.names == pp[2]) %>%
        select(ensg) %>% unique()

    n <- length(intersect(g1$ensg, g2$ensg))
    m <- length(union(g1$ensg, g2$ensg))

    data.frame(n, m, t1 = pp[1], t2 = pp[2])
}

tis.jaccard <- do.call(rbind, apply(tis.pairs, 2, take.jaccard)) %>%
    mutate(ti.1 = match(t1, tis.names),
           ti.2 = match(t2, tis.names))

n.tis <- length(tis.names)
J <- matrix(NA, n.tis, n.tis)
J[tis.jaccard %>% select(ti.1, ti.2) %>% as.matrix] <- tis.jaccard$n / tis.jaccard$m
J[tis.jaccard %>% select(ti.2, ti.1) %>% as.matrix] <- tis.jaccard$n / tis.jaccard$m

hh <- hclust(as.dist(1 - J))
oo <- order.optimal(as.dist(1 - J), hh$merge)

tis.names.sorted <- tis.names[oo$order]

j.df.1 <- tis.jaccard %>% select(n, m, t1, t2) %>% as.data.frame()
j.df.2 <- tis.jaccard %>% select(n, m, t2, t1) %>% as.data.frame()
colnames(j.df.1) <- c('n', 'm', 't1', 't2')
colnames(j.df.2) <- c('n', 'm', 't1', 't2')
j.df <- rbind(j.df.1, j.df.2)

j.df$t1 <- factor(j.df$t1, tis.names.sorted)
j.df$t2 <- factor(j.df$t2, tis.names.sorted)

j.diag <- data.frame(tt = tis.names)
j.diag$tt <- factor(j.diag$tt, tis.names.sorted)

p3 <- ggplot() + theme_bw() +
    geom_point(data = j.diag, aes(x = tt, y = tt, fill = tt),
               pch = 22, size = 3.5, show.legend = FALSE) +
    geom_point(data = j.df %>% filter(as.integer(t1) < as.integer(t2)),
               aes(x = t1, y = t2, size = n/m), color = 'gray20', pch = 15)

p3 <- p3 + geom_text(data = j.df %>% filter(n / m > .3, as.integer(t1) < as.integer(t2)),
               aes(x = t1, y = t2, label = round(10 * n/m)), size = 1.25, color = 'white')

my_trans <- trans_new(name = 'my', transform = function(x) round(x^2.5, 3),
                      inverse = function(x) round((x)^(1/2.5), 3))

dl.p3 <- list(dl.trans(x = x + .2), 'last.points', cex = .6, rot = 0,
              hjust = 0, vjust = 0.5)

p3 <- p3 + geom_dl(data = j.diag, aes(x = tt, y = tt, label = tt),
             method = dl.p3) +
    scale_fill_manual(values = tis.colors[oo$order, 2]) +
    scale_size_continuous('Jaccard', range = c(0, 3.5), trans = my_trans) +
    theme(axis.text.x = element_blank(), axis.title = element_blank(),
          axis.text.y = element_text(size = 8),
          legend.position = c(1, 0), legend.justification = c(1,0),
          legend.background = element_blank())

p3.1.df <- gene.count.tab %>%
    mutate(pr.genes = num.genes / num.genes.tot)

p3.1.df$tis.names <- factor(p3.1.df$tis.names, tis.names.sorted)

p3.1 <- ggplot(p3.1.df, aes(x = tis.names, y = pr.genes)) + theme_classic() +
    geom_bar(stat = 'identity', fill = 'gray70') +
    theme(axis.text.x = element_blank(), axis.title = element_blank())

pdf(file = tis.pair.file, width = 7, height = 6, useDingbats = FALSE)
grid.vcat(list(p3.1, p3), heights = c(.5,2.5))
dev.off()
