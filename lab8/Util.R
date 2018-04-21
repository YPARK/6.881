require(grid)
require(gridExtra)
require(gtable)
require(ggplot2)

match.widths.grob <- function(g.list) {

    max.width <- g.list[[1]]$widths[2:7]

    for(j in 2:length(g.list)) {
        max.width <- grid::unit.pmax(max.width, g.list[[j]]$widths[2:7])
    }

    for(j in 1:length(g.list)) {
        g.list[[j]]$widths[2:7] <- as.list(max.width)
    }
    return(g.list)
}

match.widths <- function(p.list) {
    g.list <- lapply(p.list, ggplotGrob)
    return(match.widths.grob(g.list))
}

grid.vcat <- function(p.list, ...) {
    g.list <- match.widths(p.list)
    ret <- grid.arrange(grobs = g.list, ncol = 1, newpage = FALSE, ...)
    return(ret)
}

match.heights.grob <- function(g.list, stretch = TRUE)  {
    max.height <- g.list[[1]]$heights[2:7]

    if(stretch) {
        for(j in 2:length(g.list)) {
            max.height <- grid::unit.pmax(max.height, g.list[[j]]$heights[2:7])
        }
    }

    for(j in 1:length(g.list)) {
        g.list[[j]]$heights[2:7] <- as.list(max.height)
    }

    return(g.list)
}

match.heights <- function(p.list, stretch = FALSE) {
    g.list <- lapply(p.list, ggplotGrob)
    return(match.heights.grob(g.list, stretch))
}

grid.hcat <- function(p.list, ...) {
    g.list <- match.heights(p.list, stretch = TRUE)
    ret <- grid.arrange(grobs = g.list, nrow = 1, newpage = FALSE, ...)
    return(ret)
}

