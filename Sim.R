.rnorm <- function(n1, n2) {
    return(matrix(rnorm(n1 * n2), nrow = n1, ncol = n2))
}

sample.tissues <- function(X,
                           n.tis,
                           n.causal.snp,
                           n.causal.tis,
                           K,
                           h2,
                           rseed) {
    
    n <- dim(X)[1]
    p <- dim(X)[2]

    set.seed(rseed)

    causal.snp <- NULL
    causal.pairs <- NULL

    eta <- matrix(0, nrow = n, ncol = n.tis)

    causal.tis <- sample(n.tis, n.causal.tis)
    tis.assign <- sample(K, n.causal.tis, replace = TRUE)

    for(k in unique(tis.assign)) {

        .snp <- sample(p, n.causal.snp)

        .tis <- causal.tis[which(tis.assign == k)]

        if(length(.tis) == 0){ next }

        .pairs <- data.frame(expand.grid(.snp, .tis), k)

        .X <- X[, .snp, drop=FALSE]

        .left <- .rnorm(n.causal.snp, 1)
        .right <- .rnorm(length(.tis), 1)

        .effect <- matrix(.left %*% t(.right),
                          nrow = dim(.left)[1],
                          ncol = dim(.right)[1])

        if(n.causal.snp > 1) {
            .effect <- as.matrix(scale(.effect))
        }

        eta[, .tis] <- eta[, .tis] + .X %*% .effect

        causal.snp <- c(causal.snp, .snp)
        causal.pairs <- rbind(causal.pairs, .pairs)
    }

    causal.snp <- unique(causal.snp)
    causal.pairs <- unique(causal.pairs)

    t.noise <- n.tis - length(causal.tis)

    if(t.noise > 0) {
        s0 <- mean(apply(eta[, causal.tis, drop = FALSE], 2, sd))
        eta[, -causal.tis] <- .rnorm(n, t.noise) * s0
    }

    ## 3. add noise
    v.eta <- apply(eta, 2, var, na.rm = TRUE)
    v.err <- v.eta * (1/h2 - 1)
    sd.err <- sqrt(v.err)

    Y <- eta + sweep(.rnorm(n, n.tis), 2, sd.err, `*`)

    colnames(causal.pairs) <- c('snp', 'tis', 'rank')
    list(output = Y, causal.pairs = causal.pairs, causal.tis = causal.tis, eta = eta)
}

add.missing.values <- function(Y, matched.tab) {
    require(dplyr)

    Y.obs <- matrix(NA, nrow = dim(Y)[1], ncol = dim(Y)[2])
    .idx <- as.matrix(matched.tab %>% select(ind.pos, tis.pos))
    Y.obs[.idx] <- Y[.idx]

    return(Y.obs)
}

get.pve <- function(y, beta, x){
    var(x %*% beta, na.rm=T) / var(y, na.rm=T)
}

get.pve.no.y <- function(beta, x){
    var(x %*% beta, na.rm=T)
}

get.lm.pve <- function(y, xx){
    .anova <- anova(lm(y~ ., data = data.frame(y, x = xx)))
    .ssq <- .anova$"Sum Sq"
    .p <- length(.ssq)
    sum(.ssq[-.p]) / sum(.ssq)
}
