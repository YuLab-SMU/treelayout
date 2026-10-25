#' Tree layout names
#'
#' Single source of truth for layout names and layout groupings.
#' Kept as a plain list rather than a registry: the groupings below already
#' exist scattered across 'ggtree'; this collects them in one place.
#'
#' @keywords internal
.layout_names <- list(
    ## every layout `tree_layout()` accepts by name
    all = c("rectangular", "slanted", "fan", "circular",
            "inward_circular", "radial", "unrooted", "equal_angle",
            "daylight", "tree_and_leaf", "dendrogram", "ape",
            "ellipse", "roundrect", "tidy"),

    ## layouts computed by `layout.unrooted()`; previously written out in
    ## full at a dozen places in ggtree
    unrooted = c("equal_angle", "daylight", "ape", "tree_and_leaf"),

    ## layout name used by tanglegram, not reachable via `ggtree(layout=)`
    internal = c("tanglegram")
)


#' Compute node coordinates for a tree
#'
#' Turn a tree into a table of node coordinates, one row per node. This is the
#' layout engine behind 'ggtree', extracted so the coordinates can be computed
#' without any plotting dependency.
#'
#' @param model a tree object: `phylo`, `treedata`, `multiPhylo`, or anything
#'   with an `as.phylo()` method.
#' @param data unused; kept so the signature matches `ggplot2::fortify()`,
#'   which callers such as 'ggtree' pass through.
#' @param layout layout name or a function
#'   taking an `igraph` object and returning an n-by-2 matrix of coordinates.
#' @param ladderize logical, ladderize the tree before layout.
#' @param right logical, ladderize to the right.
#' @param branch.length variable used to scale branch lengths; `"none"` gives a
#'   cladogram.
#' @param mrsd most recent sampling date, for time-scaled trees.
#' @param as.Date logical, whether `mrsd` should be treated as a Date.
#' @param yscale variable used to scale tips along y; `"none"` for standard.
#' @param root.position x position of the root.
#' @param ... passed on to the layout function.
#'
#' @return A `tbl_df` with one row per node and columns `node`, `parent`,
#'   `x`, `y`, `branch`, `angle`, `isTip` and `label`, plus any columns carried
#'   over from `treedata`. For a `multiPhylo` input the results are stacked and
#'   an `.id` column identifies the tree. The attribute `"layout"` records the
#'   layout used (`"external"` when `layout` was a function).
#'
#' @examples
#' tr <- ape::rtree(10)
#' head(tree_layout(tr))
#' head(tree_layout(tr, layout = "circular"))
#'
#' @importFrom ape ladderize
#' @importFrom dplyr full_join mutate
#' @importFrom tibble tibble
#' @importFrom tidytree as_tibble as.phylo Nnode
#' @export
tree_layout <- function(model, data = NULL,
                        layout        = "rectangular",
                        ladderize     = TRUE,
                        right         = FALSE,
                        branch.length = "branch.length",
                        mrsd          = NULL,
                        as.Date       = FALSE,
                        yscale        = "none",
                        root.position = 0,
                        ...) {

    ## a list of trees: lay out each and stack the results
    if (inherits(model, c("multiPhylo", "treedataList"))) {
        return(.layout_multi(model, data, layout, ladderize, right,
                             branch.length, mrsd, as.Date, yscale,
                             root.position, ...))
    }

    ## a function means an external layout (e.g. igraph::layout_with_kk)
    external <- NULL
    if (is.function(layout)) {
        external <- layout
        layout <- "equal_angle"
    } else {
        layout <- match.arg(layout, layout_names())
    }

    x <- as.phylo(model)
    if (ladderize == TRUE) {
        x <- ladderize(x, right = right)
    }

    if (!is.null(x$edge.length) && anyNA(x$edge.length)) {
        warning("'edge.length' contains NA values...\n",
                "## setting 'edge.length' to NULL automatically when plotting the tree...")
        x$edge.length <- NULL
    }

    if (!is.null(external)) {
        res <- .layout_external(model, x, external, root.position, ...)
    } else if (layout %in% .layout_names$unrooted) {
        res <- layout.unrooted(model, layout.method = layout,
                               branch.length = branch.length, ...)
    } else {
        if (layout != "tidy")
            ypos <- getYcoord(x)
        N <- Nnode(x, internal.only = FALSE)
        if (is.null(x$edge.length) || branch.length == "none") {
            if (layout == "slanted") {
                sbp <- .convert_tips2ancestors_sbp(x, include.root = TRUE)
                xpos <- getXcoord_no_length_slanted(sbp)
                ypos <- getYcoord_no_length_slanted(sbp)
            } else {
                xpos <- getXcoord_no_length(x)
            }
        } else {
            xpos <- getXcoord(x)
        }

        if (layout == "tidy")
            ypos <- getYcoord_tidy(x, xpos)

        xypos <- tibble(node = 1:N, x = xpos + root.position, y = ypos)

        df <- as_tibble(model) |>
            mutate(isTip = !.data$node %in% .data$parent)

        res <- full_join(df, xypos, by = "node")
    }

    res <- calculate_branch_mid(res, layout = layout)

    if (!is.null(mrsd)) {
        res <- scaleX_by_time_from_mrsd(res, mrsd, as.Date)
    }

    if (layout == "slanted") {
        res <- add_angle_slanted(res)
    } else {
        res <- calculate_angle(res)
    }
    res <- scaleY(as.phylo(model), res, yscale, layout, ...)
    res <- adjust_hclust_tip_edge_len(res, x, layout, branch.length)
    class(res) <- c("tbl_tree", class(res))
    attr(res, "layout") <- if (is.null(external)) layout else "external"
    return(res)
}


#' Coordinates from an external layout function
#'
#' Runs `fn` on the tree converted to an `igraph` object and returns a table
#' ready to join with the usual node columns.
#'
#' @param model the original tree, used for the node table.
#' @param x the `phylo` version, already ladderized.
#' @param fn the layout function.
#' @param root.position x offset applied to the result.
#' @param ... passed to `fn`.
#' @return a `tbl_df` with `node`, `x`, `y`.
#' @noRd
.layout_external <- function(model, x, fn, root.position, ...) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Layout functions need the 'igraph' package.\n",
             "-> install it with install.packages(\"igraph\")", call. = FALSE)
    }
    g <- ape::as.igraph.phylo(x, use.labels = FALSE)
    xy <- fn(g, ...)
    xy <- as.matrix(xy)[, 1:2, drop = FALSE]
    if (nrow(xy) != (Nnode(x, internal.only = FALSE))) {
        stop("Layout function returned ", nrow(xy), " coordinates for a tree with ",
             Nnode(x, internal.only = FALSE), " nodes.", call. = FALSE)
    }
    df <- as_tibble(model) |>
        mutate(isTip = !.data$node %in% .data$parent)
    xypos <- tibble(node = seq_len(nrow(xy)),
                    x = xy[, 1] + root.position,
                    y = xy[, 2])
    full_join(df, xypos, by = "node")
}


#' Lay out each tree in a list and stack the results
#'
#' @param ... arguments forwarded to [tree_layout()].
#' @return a `tbl_df` with an extra `.id` column identifying the tree.
#' @noRd
.layout_multi <- function(model, data, layout, ladderize, right,
                          branch.length, mrsd, as.Date, yscale,
                          root.position, ...) {
    df.list <- lapply(model, function(x) {
        tree_layout(x, data, layout = layout, ladderize = ladderize,
                    right = right, branch.length = branch.length,
                    mrsd = mrsd, as.Date = as.Date, yscale = yscale,
                    root.position = root.position, ...)
    })
    ids <- if (is.null(names(model))) paste0("Tree #", seq_along(model)) else names(model)
    names(df.list) <- ids

    ## rbind() keeps each part's class, which stacks "tbl_tree" once per tree
    df <- do.call(rbind, df.list)
    class(df) <- c("tbl_tree", "tbl_df", "tbl", "data.frame")
    df$.id <- factor(rep(ids, times = vapply(df.list, nrow, integer(1))),
                     levels = ids)
    ## rbind() kept the layout attribute of the first tree only
    attr(df, "layout") <- if (is.function(layout)) "external" else layout
    df
}


#' Layout names accepted by [tree_layout()]
#'
#' @return a character vector of layout names.
#' @examples
#' layout_names()
#' @export
layout_names <- function() {
    .layout_names$all
}


#' Read the layout recorded on a layout result
#'
#' @param x a table returned by [tree_layout()], or a plot object whose `$data`
#'   carries the layout attribute.
#' @return the layout name, or `NULL`.
#' @examples
#' tr <- ape::rtree(10)
#' get_layout(tree_layout(tr, layout = "circular"))
#' @export
get_layout <- function(x) {
    if (is.data.frame(x)) return(attr(x, "layout"))
    attr(x$data, "layout")
}
