## Frozen copies of the slanted coordinate path that .slanted_coords() replaced.
##
## They build an n_tip x n_node indicator matrix and then summarise it, which is
## the definition .slanted_coords() must reproduce -- but it is O(n^2) in
## memory, so it lives here rather than in the package. Test-only.

.internal_ancestor_ref <- function(.data, .node, all.nodes, type = "all",
                                   include.self = TRUE) {
    .internal_anc <- switch(type, all = tidytree::ancestor, parent = tidytree::parent)
    x <- .internal_anc(.data = .data, .node = .node)
    if (include.self) {
        x <- c(x, .node)
    }
    x <- all.nodes %in% x
    x
}

.convert_tips2ancestors_sbp_ref <- function(tree, include.root = FALSE,
                                            type = "all", include.self = TRUE) {
    all.nodes <- seq_len(max(tree$edge))
    if (!include.root) {
        all.nodes <- setdiff(all.nodes, ape::Ntip(tree) + 1L)
    }
    tip.nodes <- seq_len(ape::Ntip(tree))
    sbp <- lapply(tip.nodes, .internal_ancestor_ref, .data = tree,
                  all.nodes = all.nodes, type = type, include.self = include.self)
    sbp <- stats::setNames(sbp, tip.nodes)
    sbp <- do.call(rbind, sbp)
    colnames(sbp) <- all.nodes
    sbp
}

getXcoord_no_length_slanted_ref <- function(x) {
    x <- -colSums(x)
    x <- unname(x[order(as.numeric(names(x)))])
    x + max(abs(x))
}

getYcoord_no_length_slanted_ref <- function(y) {
    y <- seq_len(nrow(y)) * y
    y[y == 0] <- NA
    y <- colMeans(y, na.rm = TRUE)
    unname(y[order(as.numeric(names(y)))])
}
