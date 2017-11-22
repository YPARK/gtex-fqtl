source('Util.R')
require(Matrix)
require(methods)

################################################################
## summary-based imputation
## 
## numerator :
##
##   z.GWAS (marginal)' * z.QTL (polygenic)
##
## denominator:
##
##   eta.QTL = D * V' * z.QTL
##   eta.GWAS = inv(D) * V' * z.GWAS
##
##   1. eta.QTL' * eta.QTL
##   2. eta.GWAS' * eta.GWAS   (These were too strong)
##   3. # SNPs                 (These were too strong)
##

func.NWAS <- function(.qtl.z, gwas.z, V.t, D) {
    .num <- sum(.qtl.z * gwas.z)
    .eta.QTL <- sweep(sweep(V.t, 2, .qtl.z, `*`), 1, D, `*`)
    ## .eta.GWAS <- sweep(sweep(V.t, 2, gwas.z, `*`), 1, D, `/`)
    .denom <- sum(.eta.QTL^2)
    ## + sum(.eta.GWAS^2) + length(.qtl.z)
    .sd <- sqrt(.denom)
    return(data.frame(z = .num/.sd, theta = .num, theta.se = .sd))
}

## Calculate SVD of reference genotype matrix X
func.LD.svd <- function(X, eig.tol = 1e-8, normalize = FALSE){
    ## use centered covariance matrix
    if(normalize){
        x.safe <- scale(X, scale = TRUE, center = TRUE)
    } else {
        x.safe <- scale(X, scale = FALSE, center = TRUE)
    }
    n.ind <- apply(!is.na(x.safe), 2, sum)
    x.safe <- sweep(x.safe, 2, sqrt(pmax(n.ind, 1)), `/`)
    x.safe[is.na(x.safe)] <- 0

    ## prevent underflow of Lapack dgesdd, (-5, 5) seems safe
    n.max <- 5 / (max(abs(x.safe)) + 1e-8)
    svd.out <- svd(x.safe * n.max)
    .valid <- which(svd.out$d^2 > eig.tol)

    ## scale back eigen values
    d <- svd.out$d[.valid] / n.max
    V.t <- t(svd.out$v[, .valid, drop = FALSE])

    log.msg('SVD finished: %d x %d V.t', dim(V.t)[1], dim(V.t)[2])

    return(list(V.t = V.t, d = d))
}

## Correct GWAS z-score bias
func.correct.bias <- function(z, V.t, D){
    V.t1 <- matrix(apply(V.t, 1, sum), ncol=1)
    num <- sweep(t(V.t), 2, D^2, `*`) %*% V.t1 * sum(z)
    denom <- sum((D * V.t1)^2)
    return(z - num/denom)
}

func.blk.ind <- function(n.qtl, n.blk) {

    ret <- sparseMatrix(i = as.integer(1:(n.qtl*n.blk)),
                        j = as.integer(as.vector(sapply(1:n.blk, function(x) rep(x, n.qtl)))),
                        x = 1)
    
    log.msg('Made %d x %d block index matrix\n', n.qtl*n.blk, n.blk)

    return(ret)
}

## faster QTL permutation
func.NWAS.qtl.perm <- function(qtl.z, gwas.z, V.t, D, blk.ind, DEBUG = FALSE) {

    n.blk <- dim(blk.ind)[2]
    n.tot.snp <- length(gwas.z)
    n.qtl <- length(qtl.z)
    qtl.idx <- sample(n.tot.snp, n.blk * n.qtl, replace = TRUE)
    
    stat.num <- sweep(matrix(gwas.z[qtl.idx], nrow = n.qtl, byrow = FALSE), 1, qtl.z, `*`)
    stat.num <- apply(stat.num, 2, sum)
    
    eta.QTL.perm <- 
        sweep(sweep(V.t[, qtl.idx, drop = FALSE], 1, D, `*`),
              2, rep(qtl.z, n.blk), `*`)
    
    stat.denom <- matrix(apply(eta.QTL.perm^2, 2, sum), nrow = 1) %*% blk.ind
    
    stat.num <- as.numeric(stat.num)
    stat.denom <- as.numeric(stat.denom)

    if(DEBUG){
        ## confirmed: slow permutation by permutation statistic
        .test <- matrix(qtl.idx, nrow = n.qtl, byrow = FALSE)
        .func <- function(.idx) {
            func.NWAS(qtl.z, gwas.z[.idx], V.t[, .idx, drop = FALSE], D)
        }
        test.temp <- do.call(rbind, apply(.test, 2, .func))
        cat(sum(abs(stat.num/sqrt(stat.denom) - test.temp[, 1])))
    }

    return(data.frame(z = stat.num/sqrt(stat.denom),
                      theta = stat.num,
                      theta.se = sqrt(stat.denom)))
}
