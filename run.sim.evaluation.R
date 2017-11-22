#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 2) {
    q()
}

sim.tag <- argv[1] # e.g., '19077/1/0.1/1/10/1001'
out.file <- argv[2]

data.name <- function(ext) {
    paste('scratch/simulation/data/', sim.tag, '-temp', ext, sep = '')
}

result.name <- function(method, ext) {
    paste('scratch/simulation/result/', method, '/', sim.tag, '-', ext, sep = '')
}

library(dplyr)
library(PRROC)
library(reshape2)
source('Util.R')

get.power <- function(score, lab, fdr.cutoff = 0.01) {
    pr.out <- pr.curve(score[lab == 1], score[lab == 0], curve = TRUE)
    .idx <- pr.out$curve[, 2] >= (1 - fdr.cutoff)
    if(sum(.idx) == 0) return(0)
    .temp <- pr.out$curve[.idx, 1]
    return(.temp[1])
}

get.auprc <- function(score, lab) {
    pr.out <- pr.curve(score[lab == 1], score[lab == 0])
    return(pr.out$auc.davis.goadrich)
}

options(stringsAsFactors = FALSE)

causal.tab <- read.table(data.name('.causal.txt.gz'))

#####################################
## comparison of tissue prediction ##
#####################################

phenix.pve <- read.table(result.name('phenix', 'phenix.pve.gz')) %>%
    rename(tis = V1, score = V2) %>%
        mutate(method = 'phenix')

bslmm.pve <- read.table(result.name('bslmm', 'bslmm.hyper.gz')) %>%
    filter(V1 == 'pve') %>%
        select(V8, V2) %>% rename(tis = V8, score = V2) %>%
            mutate(method = 'bslmm-pve')

bslmm.pge <- read.table(result.name('bslmm', 'bslmm.hyper.gz')) %>%
    filter(V1 == 'pge') %>%
        select(V8, V2) %>% rename(tis = V8, score = V2) %>%
            mutate(method = 'bslmm-pge')

fqtl.pve <- read.table(result.name('fqtl', 'fqtl.pve.gz')) %>%
    mutate(method = paste(V4, V5, sep = '-')) %>%
        rename(tis = V1, score = V2) %>% select(tis, score, method)

zfqtl.pve <- read.table(result.name('zfqtl', 'zfqtl.pve.gz')) %>%
    mutate(method = paste(V4, V5, sep = '-')) %>%
        rename(tis = V1, score = V2) %>% select(tis, score, method)

sqtl.pve <- read.table(result.name('sqtl', 'sqtl.pve.gz')) %>%
    rename(tis = V1, score = V2, method = V3)

metasoft.m.tis <- read.table(result.name('metasoft', 'metasoft.tis.gz')) %>%
    rename(tis = V2, score = V4) %>% mutate(method = 'metasoft-m') %>%
        select(tis, score, method)

metasoft.p.tis <- read.table(result.name('metasoft', 'metasoft.tis.gz')) %>%
    rename(tis = V2) %>% mutate(score = -log10(1e-10 + V3), method = 'metasoft-p') %>%
        select(tis, score, method)

mashr.tis <- read.table(result.name('mashr', 'mashr.tis.gz')) %>%
    rename(tis = V2) %>% mutate(score = -log10(1e-10 + V4), method = 'mashr') %>%
        select(tis, score, method)

glmnet.pve <- read.table(result.name('glmnet', 'glmnet.pve.gz')) %>%
    mutate(method = paste(V3, V4, sep = '-')) %>%
        rename(tis = V1, score = V2) %>% select(tis, score, method)

combined <- rbind(fqtl.pve, zfqtl.pve, sqtl.pve, phenix.pve, bslmm.pve, bslmm.pge,
                  metasoft.m.tis, metasoft.p.tis, mashr.tis, glmnet.pve) %>%
                      mutate(lab = ifelse(tis %in% causal.tab[, 2], 1, 0))

tis.out <- combined %>% na.omit() %>% group_by(method) %>%
    summarize(tis.auprc = get.auprc(score, lab),
              tis.power.01 = get.power(score, lab, fdr.cutoff = 0.01),
              tis.power.10 = get.power(score, lab, fdr.cutoff = 0.10))


##################################
## comparison of snp prediction ##
##################################

phenix.snp <- read.table(result.name('phenix', 'phenix.snp.gz')) %>%
    rename(snp = V1) %>%
        mutate(score = abs(V3), method = paste('phenix', V4, sep = '-')) %>%
            select(snp, score, method)

glmnet.snp <- read.table(result.name('glmnet', 'glmnet.snp.gz')) %>%
    rename(snp = V1) %>%
        mutate(score = abs(V3), method = paste('glmnet', V4, sep = '-')) %>%
            select(snp, score, method)

bim.tab <- read.table(data.name('.bim'))

fqtl.snp <- read.table(result.name('fqtl', 'fqtl.snp.gz')) %>%
    mutate(method = paste('fqtl', V6, sep = '-'), snp = bim.tab[V1,2],
           score = abs(V4)) %>%
               select(snp, score, method)

zfqtl.snp <- read.table(result.name('zfqtl', 'zfqtl.snp.gz')) %>%
    mutate(method = paste('zfqtl', V6, sep = '-'), snp = bim.tab[V1,2],
           score = abs(V4)) %>%
               select(snp, score, method)

sqtl.snp <- read.table(result.name('sqtl', 'sqtl.snp.gz'))
colnames(sqtl.snp) <- c('snp', 1:44, 'method')
sqtl.snp <- melt(sqtl.snp, id.vars = c('snp', 'method'),
                 variable.name = 'tissue',
                 value.name = 'beta')

sqtl.snp <- sqtl.snp %>% mutate(score = abs(beta)) %>%
    group_by(snp, method) %>% slice(which.max(score)) %>%
        select(snp, score, method)

bslmm.snp.full <- read.table(result.name('bslmm', 'bslmm.snp.gz'), header = TRUE) %>%
    mutate(score = abs(dense + sparse * pip), method = 'bslmm.full') %>%
        rename(snp = rs) %>% group_by(snp, method) %>% slice(which.max(score)) %>%
            select(snp, score, method)

bslmm.snp.sparse <- read.table(result.name('bslmm', 'bslmm.snp.gz'), header = TRUE) %>%
    mutate(score = abs(sparse * pip), method = 'bslmm.sparse') %>%
        rename(snp = rs) %>% group_by(snp, method) %>% slice(which.max(score)) %>%
            select(snp, score, method)

metasoft.snp <- read.table(result.name('metasoft', 'metasoft.snp.gz')) %>%
    mutate(score = -log10(1e-8 + V2), method = 'metasoft') %>% rename(snp = V1) %>%
        select(snp, score, method)

mashr.snp.lfsr <- read.table(result.name('mashr', 'mashr.snp.gz')) %>%
    mutate(score = -log10(1e-8 + V3), method = 'mashr-lfsr') %>% rename(snp = V1) %>%
        select(snp, score, method)

mashr.snp.lfdr <- read.table(result.name('mashr', 'mashr.snp.gz')) %>%
    mutate(score = -log10(1e-8 + V2), method = 'mashr-lfdr') %>% rename(snp = V1) %>%
        select(snp, score, method)

combined <- list(fqtl.snp, zfqtl.snp, sqtl.snp, phenix.snp, bslmm.snp.full, bslmm.snp.sparse,
                 metasoft.snp, mashr.snp.lfsr, mashr.snp.lfdr, glmnet.snp)

combined <- do.call(rbind, lapply(combined, as.data.frame)) %>%
    mutate(lab = ifelse(snp %in% causal.tab[, 4], 1, 0))

snp.out <- combined %>% na.omit() %>% group_by(method) %>%
    summarize(tis.auprc = get.auprc(score, lab),
              tis.power.01 = get.power(score, lab, fdr.cutoff = 0.01),
              tis.power.10 = get.power(score, lab, fdr.cutoff = 0.10))

tot.out <- as.data.frame(rbind(tis.out %>% mutate(task = 'tissue'),
                               snp.out %>% mutate(task = 'snp')))

tot.out <- data.frame(tot.out, tag = gsub('/', '\t', sim.tag))

write.tsv(tot.out, file = gzfile(out.file))
