#' Calculate angle coordinates for circular layouts
#'
#' Adds an `angle` column holding the label angle for every node, derived from
#' the node's position and its parent.
#'
#' @param data a node table with `x`, `y`, `node` and `parent` columns.
#' @return the node table with an `angle` column added.
#' @examples
#' d <- tree_layout(ape::rtree(6), layout = "circular")
#' calculate_angle(d)
#' @export
calculate_angle <- function(data) {
    if (!"y" %in% names(data)) {
        stop("Input data must contain 'y' column")
    }
    y_range <- diff(range(data$y, na.rm = TRUE))
    if (y_range == 0) {
        warning("Y range is 0, setting all angles to 0")
        data$angle <- 0
    } else {
        y_max <- max(data$y, na.rm = TRUE)
        y_min <- min(data$y, na.rm = TRUE)
        if (y_max > y_range + 1) {
            ## The formula below assumes the largest y is y_range + 1, which
            ## holds while y runs 1..n -- what every layout here produces, and
            ## also what the unrooted layouts give when their y is negative
            ## (their y_max stays below y_range + 1).
            ##
            ## ggtree's scaleClade() rescales one clade and shifts the rest,
            ## lifting the minimum above 1. y_max then exceeds y_range + 1 and
            ## the angle passed 360 (473 at 20 tips), so the labels on that side
            ## pointed the wrong way. #647
            ##
            ## Shifting y onto 1..n only in that case keeps every layout's
            ## angles exactly as they were.
            data$angle <- 360 * (data$y - y_min + 1) / (y_range + 1)
        } else {
            data$angle <- 360 / (y_range + 1) * data$y
        }
    }
    return(data)
}


#' Scale Y coordinates based on specified variable
#'
#' @param phylo Phylogenetic tree object
#' @param df Data frame with tree data
#' @param yscale Variable name to use for scaling
#' @param layout Tree layout
#' @param ... Additional arguments passed to scaling functions
#' @return Modified data frame with scaled y coordinates
#' @noRd
scaleY <- function(phylo, df, yscale, layout, ...) {
    if (yscale == "none") {
        return(df)
    }
    if (!yscale %in% colnames(df)) {
        warning("yscale variable '", yscale, "' is not available in data frame")
        return(df)
    }

    if (is.numeric(df[[yscale]])) {
        y <- getYcoord_scale_numeric(phylo, df, yscale, ...)
    } else {
        y <- getYcoord_scale_category(phylo, df, yscale, ...)
    }

    df[, "y"] <- y
    return(df)
}


#' Adjust tip edge lengths for hclust objects
#'
#' @param df Data frame with tree data
#' @param phylo Phylogenetic tree object
#' @param layout Tree layout
#' @param branch.length Branch length specification
#' @return Modified data frame with adjusted tip edge lengths
#' @noRd
adjust_hclust_tip_edge_len <- function(
    df,
    phylo,
    layout,
    branch.length = "branch.length"
) {
    # Determine tip.edge.len based on object type
    if (inherits(phylo, 'treedata')) {
        tip.edge.len <- attr(phylo@phylo, 'tip.edge.len')
    } else {
        tip.edge.len <- attr(phylo, 'tip.edge.len')
    }

    # Apply adjustments if conditions met
    if (
        !is.null(tip.edge.len) &&
            branch.length != 'none' &&
            ## layouts where the tip edge length is meaningful, i.e. the
            ## ones that hang the tips below the rest of the tree
            layout %in% c('rectangular', 'dendrogram', 'roundrect')
    ) {
        mx <- max(df$x, na.rm = TRUE)
        df$x <- df$x - mx
        df$branch <- df$branch - mx
        df[df$isTip, "x", drop = TRUE] <- tip.edge.len
        attr(df, 'revts.done') <- TRUE
    }
    return(df)
}
