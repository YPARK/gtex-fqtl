#!/usr/bin/env Rscript
## show power calculation results

options(stringsAsFactors = FALSE)

library(dplyr)
library(ggplot2)

dat.tab <- read.table('result/simulation/power.txt.gz')

colnames(dat.tab) <- c('method', 'auprc', 'power.1', 'power.10', 'task',
                       'gene', 'rank', 'h2', 'n.snp', 'n.tis', 'rseed')

dat.tab <- dat.tab %>% na.omit()

tis.n <- sort(unique(dat.tab$n.tis))
dat.tab$n.tis <- factor(dat.tab$n.tis, tis.n, paste('#tis of action =', tis.n))

###################################
## draw tissue of action results ##
###################################

str.1 <- function(s) strsplit(as.character(s), split='-')[[1]][1]
str.2 <- function(s) strsplit(as.character(s), split='-')[[1]][2]

rename.method <- function(s) {
    .s <- strsplit(as.character(s), split = '-')[[1]]
    .m <- .s[1]
    .s2 <- .s[2]

    ## alpha < 1 : elastic-net
    ## alpha == 1 : lasso
    if(.m == 'glmnet') {
        .s2 <- as.numeric(.s2)
        return(ifelse(.s2 < 1, paste('elastic net',.s2,sep='-'),
                      paste('lasso','',sep='-')))
    } else if(.m == 'fqtl') {
        .s2 <- as.numeric(.s2)
        return(ifelse(.s2 == 5, 'fqtl', NA))
    } else if(.m == 'zfqtl') {
        .s2 <- as.numeric(.s2)
        return(ifelse(.s2 == 5, 'zfqtl', NA))
    } else if(.m == 'bslmm') {
        return(ifelse(.s2 == 'pge', 'bslmm', NA))
    }
    return(s)
}

tis.power.tab <- dat.tab %>%
    filter(task == 'tissue') %>%
    group_by(method, rank, h2, n.snp, n.tis) %>%
            summarize(power.mean = mean(power.1),
                      power.se = sd(power.1)/sqrt(length(power.1))) %>%
    as.data.frame() %>%
    mutate(method = sapply(method, rename.method)) %>%
    filter(!is.na(method)) %>%
    as.data.frame() %>%
    mutate(method.setting = sapply(method, str.2),
           method = sapply(method, str.1)) %>%
    group_by(method, rank, h2, n.snp, n.tis) %>%
    slice(which.max(power.mean)) %>%
    select(-method.setting) %>%
    as.data.frame()

################################################################
## 1a. for fixed number of SNPs show different ranks
## FQTL paper Fig.2(a)

plt.power.aes <- aes(x = h2, y = power.mean, color = method,
                     lty = method, shape = method, 
                     ymin = power.mean - 2 * power.se,
                     ymax = power.mean + 2 * power.se)

plt.func <- function(.nsnp, .data, .aes, .ylab,
                      .methods, .method.names, .colors, .lty, .pch) {
    
   
    .tab <- .data %>% filter(method %in% .methods, n.snp == .nsnp)
    
    .kk <- sort(unique(.tab$rank))
    .tab$rank <- factor(.tab$rank, .kk, paste('#regulatory mech=',.kk))
    
    .tab$method <- factor(.tab$method, .methods, .method.names)
    

    plt <-
        ggplot(.tab %>% na.omit(), .aes) +
            geom_errorbar(width = 0.005, size = 0.25) +
                geom_line(size = .75) +
                    geom_point(size = 1.25, fill = 'white') +
                        facet_grid(rank ~ n.tis, scales = 'free')

    plt <- plt + theme_classic() +
        scale_color_manual(values = .colors) +
            scale_shape_manual(values = .pch) +
                scale_fill_manual(values = .colors) +
                    scale_linetype_manual(values = .lty) +
                        xlab('heritability') + ylab(.ylab)

    plt <- plt +
        theme(panel.border = element_rect(color = 'gray80', fill = NA,
                  size = 0.25),
              strip.background = element_rect(color = 'gray80', fill = 'gray80',
                  size = 0.25))

    return(plt)
}

r.dir <- 'result/simulation/figures'

n.snps <- unique(dat.tab$n.snp)
plt.list <- lapply(n.snps, plt.func,
                   .data = tis.power.tab,
                   .aes = plt.power.aes,
                   .ylab = 'power (FDR < 1%)',
                   .methods = c('fqtl', 'metasoft', 'mashr', 
                       'phenix', 'sqtl', 'bslmm',
                       'elastic net', 'lasso'),
                   .method.names = c('fQTL', 'Metasoft', 'MASH',
                       'PHENIX', 'sQTL', 'BSLMM',
                       'Elastic net', 'LASSO'),
                   .colors = c('#3333ff', 'gray40', 'gray60', '#ff1111',
                       '#6666ff', '#119911', 'gray60', 'gray40'),
                   .lty = c(1, 4, 1, 1, 2, 1, 1, 2),
                   .pch = c(21, 32, 32, 24, 32, 18, 22, 25))

plt.files <- paste(r.dir, '/fig1a-sim-tis-power-snp', n.snps, '.pdf', sep = '')
sapply(1:length(plt.files),
       function(j, ...) ggsave(plt.files[j], plot = plt.list[[j]], ...),
       width = 8, height = 5, useDingbats = FALSE)

## 1b. just compare fqtl, zQTL and phenix

plt.list <- lapply(n.snps, plt.func,
                   .data = tis.power.tab,
                   .aes = plt.power.aes,
                   .ylab = 'power (FDR < 1%)',
                   .methods = c('fqtl', 'zfqtl', 'phenix'),
                   .method.names = c('fQTL', 'fQLTL summary', 'PHENIX'),
                   .colors = c('gray40', 'gray20', '#ff3333'),
                   .lty = c(1, 4, 1),
                   .pch = c(21, 32, 32))

plt.files <- paste(r.dir, '/fig1b-sim-tis-approx-snp', n.snps, '.pdf', sep = '')
sapply(1:length(plt.files),
       function(j, ...) ggsave(plt.files[j], plot = plt.list[[j]], ...),
       width = 8, height = 5, useDingbats = FALSE)


######################
## draw snp results ##
######################

################################################################
## 2a. for a fixed number of causal SNPs show different ranks
## FQTL paper Fig.2(b)

#######################
## Power at FDR < 1% ##
#######################

snp.power.tab <- dat.tab %>%
    filter(task == 'snp') %>%
    group_by(method, rank, h2, n.snp, n.tis) %>%
            summarize(power.mean = mean(power.1),
                      power.se = sd(power.1)/sqrt(length(power.1))) %>%
    filter(method != 'bslmm.full') %>%
    as.data.frame()

snp.power.tab <- snp.power.tab %>%
    mutate(method = sapply(method, rename.method)) %>%
    filter(!is.na(method)) %>%
    as.data.frame()

snp.power.tab$method[snp.power.tab$method == 'bslmm.sparse'] <- 'bslmm'

snp.power.tab <- snp.power.tab %>%
    mutate(method.setting = sapply(method, str.2),
           method = sapply(method, str.1)) %>%
    group_by(method, rank, h2, n.snp, n.tis) %>%
    slice(which.max(power.mean)) %>%
    select(-method.setting) %>%
    as.data.frame()

plt.list <- lapply(n.snps, plt.func,
                   .data = snp.power.tab,
                   .aes = plt.power.aes,
                   .ylab = 'power (FDR < 1%)',
                   .methods = c('fqtl', 
                       'phenix', 'sqtl', 'bslmm',
                       'metasoft', 'mashr', 
                       'elastic net', 'lasso'),
                   .method.names = c('fQTL', 
                       'PHENIX', 'sQTL', 'BSLMM',
                       'Metasoft', 'MASH',
                       'Elastic net', 'LASSO'),
                   .colors = c('#3333ff', '#ff1111',
                       '#6666ff', '#119911', 'gray40', 'gray60',
                       'gray60', 'gray40'),
                   .lty = c(1, 1, 2, 1, 4, 1, 1, 2),
                   .pch = c(21, 24, 32, 18, 32, 32, 22, 25))

plt.files <- paste(r.dir, '/fig1c-sim-snp-power-snp', n.snps, '.pdf', sep = '')
sapply(1:length(plt.files),
       function(j, ...) ggsave(plt.files[j], plot = plt.list[[j]], ...),
       width = 8, height = 5, useDingbats = FALSE)


###########
## auprc ##
###########

snp.auprc.tab <- dat.tab %>%
    filter(task == 'snp') %>%
    group_by(method, rank, h2, n.snp, n.tis) %>%
            summarize(auprc.mean = mean(auprc),
                      auprc.se = sd(auprc)/sqrt(length(auprc))) %>%
    filter(method != 'bslmm.full') %>%
    as.data.frame()

snp.auprc.tab <- snp.auprc.tab %>%
    mutate(method = sapply(method, rename.method)) %>%
    filter(!is.na(method)) %>%
    as.data.frame()

snp.auprc.tab$method[snp.auprc.tab$method == 'bslmm.sparse'] <- 'bslmm'

snp.auprc.tab <- snp.auprc.tab %>%
    mutate(method.setting = sapply(method, str.2),
           method = sapply(method, str.1)) %>%
    group_by(method, rank, h2, n.snp, n.tis) %>%
    slice(which.max(auprc.mean)) %>%
    select(-method.setting) %>%
    as.data.frame()


plt.auprc.aes <- aes(x = h2, y = auprc.mean, color = method,
                     lty = method, shape = method, 
                     ymin = auprc.mean - 2 * auprc.se,
                     ymax = auprc.mean + 2 * auprc.se)

plt.list <- lapply(n.snps, plt.func,
                   .data = snp.auprc.tab,
                   .aes = plt.auprc.aes,
                   .ylab = 'AUPRC',
                   .methods = c('fqtl', 
                       'phenix', 'sqtl', 'bslmm',
                       'metasoft', 'mashr', 
                       'elastic net', 'lasso'),
                   .method.names = c('fQTL', 
                       'PHENIX', 'sQTL', 'BSLMM',
                       'Metasoft', 'MASH',
                       'Elastic net', 'LASSO'),
                   .colors = c('#3333ff', '#ff1111',
                       '#6666ff', '#119911', 'gray40', 'gray60',
                       'gray60', 'gray40'),
                   .lty = c(1, 1, 2, 1, 4, 1, 1, 2),
                   .pch = c(21, 24, 32, 18, 32, 32, 22, 25))

plt.files <- paste(r.dir, '/fig1d-sim-snp-auprc-snp', n.snps, '.pdf', sep = '')
sapply(1:length(plt.files),
       function(j, ...) ggsave(plt.files[j], plot = plt.list[[j]], ...),
       width = 8, height = 5, useDingbats = FALSE)

