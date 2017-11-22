library(MatrixEQTL)

make.s.data <- function(mat, row.names, col.names) {
    stopifnot(is.matrix(mat))
    ret <- SlicedData$new()
    ret$CreateFromMatrix(mat)
    rownames(ret) <- row.names
    colnames(ret) <- col.names
    return(ret)
}

run.matrix.qtl <- function(X, Y, C, snp.names, individuals, me.out.file) {

    n.tis <- dim(Y)[2]
    n.covar <- dim(C)[2]

    y.data <- make.s.data(t(Y), as.character(1:n.tis), individuals)
    x.data <- make.s.data(t(X), snp.names, individuals)

    if(is.null(C)){
        c.data <- SlicedData$new()
    } else {
        c.data <- make.s.data(t(C), as.character(1:n.covar), individuals)
    }

    me <- Matrix_eQTL_engine(snps = x.data,
                             gene = y.data,
                             cvrt = c.data,
                             pvOutputThreshold = 1,
                             noFDRsaveMemory = TRUE,
                             output_file_name = me.out.file)

    return(me)
}
