## Both of these were O(n^2) on ladder-like trees (depth ~ n): slanted built an
## n_tip x n_node indicator matrix, and the yscale helpers rescanned the whole
## table per node. At 2000 ladder tips slanted took 4.7 s and yscale took 4.9 s.

test_that("slanted coordinates match the indicator-matrix definition", {
    for (n in c(4, 7, 15, 40)) {
        set.seed(n)
        tr <- ape::rtree(n)
        ## the definition this replaces
        sbp <- .convert_tips2ancestors_sbp(tr, include.root = TRUE)
        expect_equal(.slanted_coords(tr)$x, getXcoord_no_length_slanted(sbp))
        expect_equal(.slanted_coords(tr)$y, getYcoord_no_length_slanted(sbp))
    }
})

test_that("slanted coordinates match on ladder and star shapes", {
    for (tr in list(ape::stree(12, type = "left"), ape::stree(9, type = "star"))) {
        sbp <- .convert_tips2ancestors_sbp(tr, include.root = TRUE)
        expect_equal(.slanted_coords(tr)$x, getXcoord_no_length_slanted(sbp))
        expect_equal(.slanted_coords(tr)$y, getYcoord_no_length_slanted(sbp))
    }
})

test_that("the slanted layout still has a sane shape", {
    tr <- ape::stree(50, type = "left")
    d <- tree_layout(tr, layout = "slanted")
    expect_equal(nrow(d), ape::Ntip(tr) + ape::Nnode(tr))
    expect_false(anyNA(d$x))
    expect_true(all(is.finite(d$y)))
    ## tips occupy one row each
    expect_equal(sort(d$y[d$isTip]), seq_len(ape::Ntip(tr)))
})

test_that("a numeric yscale fills every row", {
    ## .assign_parent_status() / .assign_child_status() used to call child() and
    ## parent() per node, each of which scanned the whole table
    tr <- ape::rtree(60)
    dat <- tibble::tibble(node = seq_len(60), val = runif(60))
    td <- new("treedata", phylo = tr, data = dat)
    d <- tree_layout(td, yscale = "val")
    expect_false(anyNA(d$y))
})

test_that("a numeric yscale fills every row on a ladder tree", {
    tr <- ape::stree(60, type = "left")
    dat <- tibble::tibble(node = seq_len(60), val = runif(60))
    td <- new("treedata", phylo = tr, data = dat)
    d <- tree_layout(td, yscale = "val")
    expect_false(anyNA(d$y))
})

test_that("children and parent indexes agree with tidytree", {
    tr <- ape::rtree(20)
    kids <- .children_index(tr)
    pars <- .parent_index(tr)
    for (v in seq_len(max(tr$edge))) {
        expect_setequal(kids[[v]], tidytree::child(as_tibble(tr), v)$node)
        p <- tidytree::parent(as_tibble(tr), v)$node
        expect_equal(pars[v], if (length(p) == 0L) 0L else p)
    }
})
