#!/usr/bin/env Rscript
## Perform confounder correction while including genotype matrix

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 4) {
    q()
}

gene.id   <- as.integer(argv[1])
do.stdize <- as.logical(argv[2])
resid.file <- argv[3]
pve.file <- argv[4]

################################################################
source('Util.R')

if(file.exists(resid.file)) {
    log.msg('File already exists')
    q()
}

if(!file.exists(glue('scratch/data/', gene.id, '/plink.bed'))) {
    log.msg('No valid SNP\n')
    q()
}


options(stringsAsFactors = FALSE)
library(fqtl)

tis.tab <- read.table(glue('scratch/data/', gene.id, '/tissues.txt.gz'))
plink <- read.plink(glue('scratch/data/', gene.id, '/plink'))
Y <- read.table(glue('scratch/data/', gene.id, '/Y.txt.gz'))

get.pve <- function(y, theta, x) {
    var(x %*% theta, na.rm = TRUE) / var(y, na.rm = TRUE)
}

write.tsv.gz <- function(mat, file.name) {
    write.tsv(signif(mat, 4), gzfile(file.name))
}

write.spike.slab <- function(effect, hdr) {
    write.tsv.gz(effect$lodds, glue(hdr, '.lodds.txt.gz'))
    write.tsv.gz(effect$theta, glue(hdr, '.theta.txt.gz'))
    write.tsv.gz(sqrt(effect$theta.var), glue(hdr, '.se.txt.gz'))
}

individuals <- read.table('data/individuals/geno.ind.txt')

## check the order of plink$FAM
fam.ind <- lapply(plink$FAM[, 1], strsplit, split = '-')
fam.ind <- sapply(fam.ind, function(x) paste(x[[1]][1:2], collapse = '-'))
stopifnot(all(fam.ind == individuals))


take.resid <- function(tis){

    vb.opt <- list(gammax = 1e3, pi.ub = -0, pi.lb = -2, vbiter = 3500, tol = 0,
                   rate = 1e-3, decay = -0.1, verbose = TRUE, out.residual = FALSE)

    tis.name <- tis.tab[tis, 2]
    cov.file <- glue('scratch/data/covariates/', tis.name, '.covariates.txt.gz')

    C <- as.matrix(read.table(cov.file))

    tis.loc <- tis.tab[tis, 1]

    log.msg('Working on %s ...\n', tis.name)

    yy <- matrix(Y[, tis.loc], ncol = 1)
    valid <- which(!is.na(yy))

    yy <- scale(yy[valid, , drop = FALSE], scale = do.stdize)
    xx <- scale(plink$BED[valid, , drop = FALSE], scale = do.stdize)
    cc <- cbind(scale(C[valid, , drop = FALSE], scale = do.stdize), 1)

    out <- fqtl.regress(y = yy, x.mean = xx, c.mean = cc, options = vb.opt)

    .pve.g <- get.pve(yy, out$mean$theta, xx)
    .pve.c <- get.pve(yy, out$mean.cov$theta, cc)

    resid <- matrix(NA, 450, 1)
    resid[valid] <- yy - cc %*% out$mean.cov$theta

    log.msg('Done -> %.2e vs %.2e, max lodds = %.2f\n\n', .pve.g, .pve.c, max(out$mean$lodds))

    return(list(resid = resid, pve.c = .pve.c, pve.g = .pve.g, sqtl.out = out))
}

n.tis <- dim(tis.tab)[1]
log.msg('# tissues = %d\n\n', n.tis)

resid.list <- lapply(1:n.tis, take.resid)

.take.pve <- function(.list) data.frame(pve.c = .list[['pve.c']], pve.g = .list[['pve.g']])
tis.names <- as.vector(sapply(tis.tab[, 2], gsub, pattern = '_Analysis', replacement = ''))

sqtl.pve.tab <- data.frame(tis.tab[, 1], tis.names, do.call(rbind, lapply(resid.list, .take.pve)))

resid.mat <- do.call(cbind, lapply(resid.list, function(.list) .list[['resid']]))
resid.mat <- scale(resid.mat, scale = do.stdize)

write.tsv.gz(resid.mat, resid.file)
write.tsv(sqtl.pve.tab, file = gzfile(pve.file))
