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
#' @param hang fraction of the tree height by which labels hang below the rest
#'   of the tree, for `dendrogram` layout of `hclust`-like input. A negative
#'   value makes them hang down from 0.
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
                        hang          = .1,
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

    ## dendrogram-like objects (dendrogram, agnes, diana, twins, linkage, hclust)
    ## have no as.phylo() method; convert them the way the layout code does
    hclust_like <- inherits(model, c("dendrogram", "linkage", "agnes", "diana", "twins"))
    if (hclust_like) {
        model <- stats::as.hclust(model)
    }
    ## `x` is the phylo used for geometry; `model` keeps any annotations, so
    ## as_tibble() is called on `model` unless it is an hclust, which has no
    ## tidytree method and whose node table comes from `x`
    if (inherits(model, "hclust")) {
        x <- .as.phylo.hclust2(model, hang = hang)
        node_table <- as_tibble(x)
    } else {
        x <- as.phylo(model)
        node_table <- as_tibble(model)
    }
    if (ladderize == TRUE) {
        x <- ladderize(x, right = right)
    }

    ## `na_free_model` is the tree the unrooted layouts get: the same one
    ## they always got, minus the NA branch lengths that would make them
    ## emit the warning a second time.
    na_free_model <- model
    if (!is.null(x$edge.length) && anyNA(x$edge.length)) {
        warning("'edge.length' contains NA values...\n",
                "## setting 'edge.length' to NULL automatically")
        x$edge.length <- NULL
        if (inherits(na_free_model, "phylo")) {
            na_free_model$edge.length <- NULL
        }
    }

    if (!is.null(external)) {
        res <- .layout_external(model, x, external, root.position, node_table, ...)
    } else if (layout %in% .layout_names$unrooted) {
        ## The unrooted layouts re-derive the tree from scratch, so they get
        ## the original `model` -- passing the ladderized `x` would change
        ## the geometry. Only the NA branch lengths are dropped first, to
        ## stop the layouts warning about the same thing again.
        res <- layout.unrooted(na_free_model, layout.method = layout,
                               branch.length = branch.length, ...)
    } else {
        if (layout != "tidy")
            ypos <- getYcoord(x)
        N <- Nnode(x, internal.only = FALSE)
        if (is.null(x$edge.length) || branch.length == "none") {
            if (layout == "slanted") {
                ## straight from the tree; the sbp indicator matrix it used to
                ## build is O(n^2) in memory
                sl <- .slanted_coords(x)
                xpos <- sl$x
                ypos <- sl$y
            } else {
                xpos <- getXcoord_no_length(x)
            }
        } else {
            xpos <- getXcoord(x)
        }

        if (layout == "tidy")
            ypos <- getYcoord_tidy(x, xpos)

        xypos <- tibble(node = 1:N, x = xpos + root.position, y = ypos)

        ## `model`, not `x`: as_tibble() on a treedata carries its annotation
        ## columns through, and the ladderized phylo would drop them
        df <- node_table |>
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
    ## The layout helpers already mark their result as a tbl_tree, so adding the
    ## class unconditionally stacks it twice ("tbl_tree" "tbl_tree" "tbl_df" ...),
    ## which breaks any exact class() comparison. Only add it when missing.
    if (!inherits(res, "tbl_tree")) {
        class(res) <- c("tbl_tree", class(res))
    }
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
#' @param node_table the node table to join the coordinates onto.
#' @param ... passed to `fn`.
#' @return the node table with `x` and `y` joined on `node`.
#' @noRd
.layout_external <- function(model, x, fn, root.position, node_table, ...) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Layout functions need the 'igraph' package.\n",
             "-> install it with install.packages(\"igraph\")", call. = FALSE)
    }
    g <- ape::as.igraph.phylo(x, use.labels = FALSE)
    xy <- fn(g, ...)
    xy <- as.matrix(xy)
    ## check the shape before subsetting: a layout function returning the wrong
    ## thing used to fail inside `[, 1:2]` with "subscript out of bounds",
    ## which says nothing about what the caller did wrong.
    if (!is.numeric(xy) || ncol(xy) < 2L) {
        stop("A layout function must return a numeric matrix with at least two ",
             "columns (x, y); got ", nrow(xy), " x ", ncol(xy), " of mode ",
             mode(xy), ".", call. = FALSE)
    }
    xy <- xy[, 1:2, drop = FALSE]
    if (nrow(xy) != (Nnode(x, internal.only = FALSE))) {
        stop("Layout function returned ", nrow(xy), " coordinates for a tree with ",
             Nnode(x, internal.only = FALSE), " nodes.", call. = FALSE)
    }
    df <- node_table |>
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
