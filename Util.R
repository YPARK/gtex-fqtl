write.tsv <- function(...) { write.table(..., col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t') }

glue <- function(...) { paste(..., sep = '') }

center <- function(...) { scale(..., center = TRUE, scale = FALSE) }

log.msg <- function(...) {
    cat('[', date() ,']', sprintf(...), file = stderr())
}

`%c%` <- function(a, b) a[, b, drop = FALSE]

`%r%` <- function(a, b) a[b, , drop = FALSE]
