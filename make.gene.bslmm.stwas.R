#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 4) q()

bslmm.weights <- argv[1] # e.g., bslmm.weights <- 'result/bslmm/1/bslmm-2-weights.gz'
data.dir <- argv[2] # e.g., data.dir <- 'scratch/data/2/'
gwas.dir <- argv[3] # e.g., gwas.dir <- 'IGAP'
out.file <- argv[4] # e.g., 'out.file.gz'

source('Util.R')

if(file.exists(out.file)){
    log.msg('File exists: %s\n\n', out.file)
    q()
}

pip.cutoff <- 0.95

n.cutoff <- 10
n.perm <- 5e7
n.blk <- 1024
n.round <- ceiling(n.perm/n.blk)

source('NWAS.R')
library(dplyr)
library(fqtl)
library(reshape2)
options(stringsAsFactors = FALSE)

plink.hdr <- glue(data.dir, '/plink')

gene.info <- read.table(glue(data.dir, '/data.gene'),
                        col.names = c('idx', 'ensg', 'chr', 'tss', 'tss.2'))
gene.info <- gene.info[2:4]
log.msg('Running on %s [%d:%d]...\n', gene.info$ensg, gene.info$chr, gene.info$tss)

## 1. read BSLMM weights
tis.info <- read.table(glue(data.dir,'/tissues.txt.gz'), col.names = c('tis', 'tis.name'))
snp.info <- read.table(glue(plink.hdr, '.bim'), col.names = c('chr', 'rs', '.', 'loc', 'qtl.a1', 'qtl.a2'))
snp.info <- snp.info %>% mutate(x.loc = 1:nrow(snp.info))

bslmm.tab <- read.table(bslmm.weights, header = TRUE)

bslmm.tab <- bslmm.tab %>% mutate(tis.name = tis.info[tis, 2]) %>%
    mutate(tis = tis.info[tis, 1]) %>%
        mutate(tis.name = gsub(tis.name, pattern = '_Analysis', replacement = ''))
    
## Use tissues with max PIP > 0.95
tissues <- bslmm.tab %>% filter(pip >= pip.cutoff) %>% select(tis) %>% unique()
tissues <- tissues[, 1]

if(length(tissues) < 1){
    system(paste('printf \"\" | gzip >', out.file))
    log.msg('No valid tissue in %s\n\n', paste(gene.info[1, ], collapse = ' '))
    q()
}

bslmm.tab <- bslmm.tab %>% filter(tis %in% tissues)

## 2. read GWAS tab
gwas.file <- glue(gwas.dir, '/chr', gene.info$chr, '.txt.gz')
plink.bim.file <- glue(plink.hdr, '.bim')
system.cmd <- paste('./make_gene_fqtl_gwas.sh', gwas.file, plink.bim.file)
gwas.txt <- system(system.cmd, intern=TRUE, ignore.stderr = TRUE)
gwas.cols <- c('chr', '.', 'loc', 'rs', 'gwas.a1', 'gwas.a2', 'gwas.z', 'gwas.theta', 'gwas.se')
gwas.tab <- read.table(text = gwas.txt, sep = '\t', col.names = gwas.cols) %>%
    mutate(chr = sapply(chr, gsub, pattern = 'chr', replacement = '')) %>%
        select(chr, loc, rs, gwas.a1, gwas.a2, gwas.z)


matched.tab <- bslmm.tab %>%
    left_join(snp.info %>% select(-chr, -loc), by = 'rs') %>%
        rename(qtl.name = rs) %>%
            left_join(gwas.tab %>% select(-chr), by = c('loc')) %>%
                na.omit() %>%
                    mutate(gwas.z.flip = if_else(qtl.a1 != gwas.a1, -gwas.z, gwas.z))


gwas.matched <- matched.tab %>% select(x.loc, gwas.z.flip) %>% unique()
gwas.z <- gwas.matched[, 'gwas.z.flip']

plink <- read.plink(plink.hdr)
LD.svd <- func.LD.svd(plink$BED, normalize = TRUE)
rm(plink)

## Select high PIP QTL
tissues <- matched.tab %>% filter(pip >= pip.cutoff) %>%
    select(tis) %>% unique()

tissues <- tissues[, 1]

if(length(tissues) < 1){
    system(paste('printf \"\" | gzip >', out.file))
    log.msg('No valid factor in %s\n\n', paste(gene.info[1, ], collapse = ' '))
    q()
}

## Permutation test
out <- NULL

for(tt in tissues) {

    obs.tab <- matched.tab %>% filter(pip >= pip.cutoff, tis == tt)
    null.tab <- matched.tab %>% filter(tis == tt)

    n.qtl <- nrow(obs.tab %>% select(qtl.name) %>% unique())
    qtl.z <- obs.tab$sparse
    blk.ind <- func.blk.ind(n.qtl, n.blk = 1024)

    vt <- LD.svd$V.t %c% obs.tab$x.loc
    dd <- LD.svd$d

    vt.null <- LD.svd$V.t %c% null.tab$x.loc
    gwas.null <- null.tab$gwas.z.flip

    obs.stat <- func.NWAS(qtl.z, obs.tab$gwas.z.flip,  vt, dd)
    z.obs.abs <- abs(obs.stat[, 'z'])

    ## adaptive permutation
    c.tot <- 0
    p.tot <- 0
    for(b in seq(1, n.round)){
        set.seed(b)

        perm.stat <- func.NWAS.qtl.perm(qtl.z, gwas.null, vt.null, dd, blk.ind)
        
        stat.z <- perm.stat[,'z']
        c.tot <- c.tot + sum(abs(stat.z) > z.obs.abs)
        p.tot <- p.tot + length(stat.z)

        if(c.tot > n.cutoff){
            break;
        }
        cat('\n', c.tot, '/', p.tot, '... ')
        cat(signif(range(stat.z), 2), '\n')
    }
    
    p.val <- (1 + c.tot)/(1 + p.tot)
    log.msg('Done : p-value = %.2e\n', p.val)

    out <- rbind(out, data.frame(obs.stat, p.val, tis = tt))
}

out <- data.frame(gene.info, gwas.dir, out) %>%
    left_join(tis.info, by = 'tis')

write.tsv(out, file = gzfile(out.file))
