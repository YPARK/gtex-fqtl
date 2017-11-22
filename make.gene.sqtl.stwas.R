#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 4) q()

sqtl.hdr <- argv[1] # e.g., sqtl.hdr <- 'result/sqtl/1/82/sqtl'
data.dir <- argv[2] # e.g., data.dir <- 'scratch/data/82/'
gwas.dir <- argv[3] # e.g., gwas.dir <- 'IGAP'
out.file <- argv[4] # e.g., 'out.file.gz'

source('Util.R')

if(file.exists(out.file)){
    log.msg('File exists: %s\n\n', out.file)
    q()
}

pip.cutoff <- 0.95
lodds.cutoff <- log(pip.cutoff) - log(1 - pip.cutoff)

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
snp.info <- read.table(glue(plink.hdr, '.bim'),
                       col.names = c('chr', 'qtl.name', '.', 'loc', 'qtl.a1', 'qtl.a2'))
snp.info <- snp.info %>% mutate(x.loc = 1:nrow(snp.info))

tis.info[, 1] <- as.integer(tis.info[, 1])
tis.info[, 2] <- sapply(tis.info[, 2], gsub, pattern = '_Analysis', replacement = '')

## linearize effects
.read.snp <- function(in.file, ...) {
    ret <- as.matrix(read.table(in.file))
    colnames(ret) <- tis.info[, 1]
    rownames(ret) <- snp.info[, 2]
    ret <- melt(ret, varnames = c('qtl.name', 'tis'), ...)
    ret$qtl.name <- as.character(ret$qtl.name)
    ret$tis <- as.character(ret$tis)
    return(ret)
}

snp.lodds <- .read.snp(glue(sqtl.hdr,'.effect.lodds.txt.gz'), value.name = 'qtl.lodds')
snp.theta <- .read.snp(glue(sqtl.hdr,'.effect.theta.txt.gz'), value.name = 'qtl.theta')
snp.se <- .read.snp(glue(sqtl.hdr,'.effect.se.txt.gz'), value.name = 'qtl.se')

sqtl.tab <- snp.theta %>% left_join(snp.se, by = c('qtl.name', 'tis')) %>%
    left_join(snp.lodds, by = c('qtl.name', 'tis'))

sqtl.tab$tis <- as.integer(sqtl.tab$tis)

## Use tissues with max PIP > 0.95
tissues <- sqtl.tab %>% filter(qtl.lodds >= lodds.cutoff) %>% select(tis) %>% unique()
tissues <- tissues[, 1]

if(length(tissues) < 1){
    system(paste('printf \"\" | gzip >', out.file))
    log.msg('No valid tissue in %s\n\n', paste(gene.info[1, ], collapse = ' '))
    q()
}

sqtl.tab <- sqtl.tab %>% filter(tis %in% tissues) %>%
    mutate(qtl.z = qtl.theta / qtl.se)

## 2. read GWAS tab
gwas.file <- glue(gwas.dir, '/chr', gene.info$chr, '.txt.gz')
plink.bim.file <- glue(plink.hdr, '.bim')
system.cmd <- paste('./make_gene_fqtl_gwas.sh', gwas.file, plink.bim.file)
gwas.txt <- system(system.cmd, intern=TRUE, ignore.stderr = TRUE)
gwas.cols <- c('chr', '.', 'loc', 'rs', 'gwas.a1', 'gwas.a2', 'gwas.z', 'gwas.theta', 'gwas.se')
gwas.tab <- read.table(text = gwas.txt, sep = '\t', col.names = gwas.cols) %>%
    mutate(chr = sapply(chr, gsub, pattern = 'chr', replacement = '')) %>%
        select(chr, loc, rs, gwas.a1, gwas.a2, gwas.z)


matched.tab <- sqtl.tab %>%
    left_join(snp.info %>% select(-chr), by = 'qtl.name') %>%
        left_join(gwas.tab %>% select(-chr), by = c('loc')) %>%
            na.omit() %>%
                mutate(gwas.z.flip = if_else(qtl.a1 != gwas.a1, -gwas.z, gwas.z))

gwas.matched <- matched.tab %>% select(x.loc, gwas.z.flip) %>% unique()
gwas.z <- gwas.matched[, 'gwas.z.flip']

plink <- read.plink(plink.hdr)
LD.svd <- func.LD.svd(plink$BED, normalize = TRUE)
rm(plink)

## Select high PIP QTL
tissues <- matched.tab %>% filter(qtl.lodds >= lodds.cutoff) %>%
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

    obs.tab <- matched.tab %>% filter(qtl.lodds >= lodds.cutoff, tis == tt)
    null.tab <- matched.tab %>% filter(tis == tt)

    n.qtl <- nrow(obs.tab %>% select(qtl.name) %>% unique())
    qtl.z <- obs.tab$qtl.z
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
