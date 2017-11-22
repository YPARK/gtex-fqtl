#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 4) q()

fqtl.hdr <- argv[1] # e.g., fqtl.hdr <- 'result/fqtl-std/51516/50/fqtl'
data.dir <- argv[2] # e.g., data.dir <- 'scratch/data/51516/'
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
library(fqtl)
library(dplyr)
library(reshape2)
options(stringsAsFactors = FALSE)

plink.hdr <- glue(data.dir, '/plink')

gene.info <- read.table(glue(data.dir, '/data.gene'),
                        col.names = c('idx', 'ensg', 'chr', 'tss', 'tss.2'))
gene.info <- gene.info[2:4]
log.msg('Running on %s [%d:%d]...\n', gene.info$ensg, gene.info$chr, gene.info$tss)

## 1. read QTL info
snp.lodds <- as.matrix(read.table(glue(fqtl.hdr,'.snp.lodds.txt.gz')))
snp.theta <- as.matrix(read.table(glue(fqtl.hdr,'.snp.theta.txt.gz')))
snp.se <- as.matrix(read.table(glue(fqtl.hdr,'.snp.se.txt.gz')))
snp.cols <- c('chr', 'qtl.name', '.', 'loc', 'qtl.a1', 'qtl.a2')
snp.info <- read.table(glue(fqtl.hdr,'.snp.info.gz'), col.names = snp.cols)
snp.info <- snp.info %>% mutate(x.loc = 1:nrow(snp.info))

tis.lodds <- read.table(glue(fqtl.hdr,'.tis.lodds.txt.gz'))
tis.theta <- as.matrix(read.table(glue(fqtl.hdr,'.tis.theta.txt.gz')))
tis.se <- as.matrix(read.table(glue(fqtl.hdr,'.tis.se.txt.gz')))
tis.info <- as.matrix(read.table(glue(fqtl.hdr,'.tis.info.gz')))

tis.info[, 1] <- as.integer(tis.info[, 1])
tis.info[, 2] <- sapply(tis.info[, 2], gsub, pattern = '_Analysis', replacement = '')

## Use factors with PIP > 0.95
factors <- which(apply(snp.lodds, 2, max) >= lodds.cutoff)

if(length(factors) < 1){
    system(paste('printf \"\" | gzip >', out.file))
    log.msg('No valid factor in %s\n\n', paste(gene.info[1, ], collapse = ' '))
    q()
}

## 2. read GWAS tab
gwas.file <- glue(gwas.dir, '/chr', gene.info$chr, '.txt.gz')
plink.bim.file <- glue(plink.hdr, '.bim')
system.cmd <- paste('./make_gene_fqtl_gwas.sh', gwas.file, plink.bim.file)
gwas.txt <- system(system.cmd, intern=TRUE, ignore.stderr = TRUE)
gwas.cols <- c('chr', '.', 'loc', 'rs', 'gwas.a1', 'gwas.a2', 'gwas.z', 'gwas.theta', 'gwas.se')
gwas.tab <- read.table(text = gwas.txt, sep = '\t', col.names = gwas.cols) %>%
    mutate(chr = sapply(chr, gsub, pattern = 'chr', replacement = '')) %>%
        select(chr, loc, rs, gwas.a1, gwas.a2, gwas.z)


qtl.z <- (snp.theta %c% factors) / (snp.se %c% factors)
colnames(qtl.z) <- factors
rownames(qtl.z) <- snp.info$qtl.name

qtl.lodds <- snp.lodds %c% factors
colnames(qtl.lodds) <- factors
rownames(qtl.lodds) <- snp.info$qtl.name

qtl.z.melt <- melt(qtl.z, varnames = c('qtl.name', 'factor'), value.name = 'qtl.z',
                   stringsAsFactors = FALSE)
qtl.z.melt$qtl.name <- as.character(qtl.z.melt$qtl.name)

qtl.lodds.melt <- melt(qtl.lodds, varnames = c('qtl.name', 'factor'), value.name = 'qtl.lodds',
                   stringsAsFactors = FALSE)

qtl.lodds.melt$qtl.name <- as.character(qtl.lodds.melt$qtl.name)

rm(qtl.z)
rm(qtl.lodds)

matched.tab <-
    qtl.z.melt %>%
        left_join(qtl.lodds.melt, by = c('qtl.name', 'factor')) %>%
            left_join(snp.info, by = 'qtl.name') %>%        
                left_join(gwas.tab %>% select(-chr), by = 'loc') %>%
                    na.omit() %>%
                        mutate(gwas.z.flip = if_else(qtl.a1 != gwas.a1, -gwas.z, gwas.z))

gwas.matched <- matched.tab %>% select(x.loc, gwas.z.flip) %>% unique()
gwas.z <- gwas.matched[, 'gwas.z.flip']

plink <- read.plink(plink.hdr)
LD.svd <- func.LD.svd(plink$BED, normalize = TRUE)
rm(plink)

## Select high PIP QTL
factors <- matched.tab %>% filter(qtl.lodds >= lodds.cutoff) %>%
    select(factor) %>% unique()

factors <- factors[, 1]

if(length(factors) < 1){
    system(paste('printf \"\" | gzip >', out.file))
    log.msg('No valid factor in %s\n\n', paste(gene.info[1, ], collapse = ' '))
    q()
}

## Permutation test
out <- NULL

for(f in factors) {

    obs.tab <- matched.tab %>% filter(qtl.lodds >= lodds.cutoff, factor == f)
    null.tab <- matched.tab %>% filter(factor == f)

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

    out <- rbind(out, data.frame(obs.stat, p.val, factor = f))
}


out <- data.frame(gene.info, gwas.dir, out)
write.tsv(out, file = gzfile(out.file))
