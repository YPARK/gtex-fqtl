
match.widths <- function(p.list) {
    require(grid)
    require(gridExtra)
    require(gtable)
    require(ggplot2)

    g.list <- lapply(p.list, ggplotGrob)
    max.width <- g.list[[1]]$widths[2:5]

    for(j in 2:length(g.list)) {
        max.width <- grid::unit.pmax(max.width, g.list[[j]]$widths[2:5])
    }

    for(j in 1:length(g.list)) {
        g.list[[j]]$widths[2:5] <- as.list(max.width)
    }
    return(g.list)
}

grid.vcat <- function(p.list, ...) {
    g.list <- match.widths(p.list)
    ret <- grid.arrange(grobs = g.list, ncol = 1, newpage = FALSE, ...)
}


match.heights <- function(p.list) {

    require(grid)
    require(gridExtra)
    require(gtable)
    require(ggplot2)

    g.list <- lapply(p.list, ggplotGrob)
    max.width <- g.list[[1]]$heights[2:5]

    for(j in 2:length(g.list)) {
        max.width <- grid::unit.pmax(max.width, g.list[[j]]$heights[2:5])
    }

    for(j in 1:length(g.list)) {
        g.list[[j]]$heights[2:5] <- as.list(max.width)
    }

    return(g.list)
}

grid.hcat <- function(p.list, ...) {
    g.list <- match.heights(p.list)
    ret <- grid.arrange(grobs = g.list, nrow = 1, newpage = FALSE, ...)
}
