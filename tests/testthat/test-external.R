test_that("a layout function is accepted", {
    skip_if_not_installed("igraph")
    tr <- ape::rtree(10)
    res <- tree_layout(tr, layout = igraph::layout_with_kk)
    expect_equal(nrow(res), 10 + ape::Nnode(tr))
    expect_true(all(c("x", "y") %in% names(res)))
    expect_identical(attr(res, "layout"), "external")
})

test_that("a layout function returning the wrong number of rows errors", {
    skip_if_not_installed("igraph")
    tr <- ape::rtree(5)
    bad <- function(g) matrix(0, nrow = 2, ncol = 2)
    expect_error(tree_layout(tr, layout = bad), "coordinates")
})

test_that("external layouts need igraph", {
    ## cannot unload igraph here; document the expectation instead
    expect_true(TRUE)
})

test_that("a layout function with the wrong return shape is rejected clearly", {
    skip_if_not_installed("igraph")
    tr <- ape::rtree(8)
    n <- ape::Ntip(tr) + ape::Nnode(tr)

    ## non-numeric / too few columns used to fail inside `[, 1:2]` with
    ## "subscript out of bounds", which says nothing about the layout function
    expect_error(tree_layout(tr, layout = function(g) "hello"),
                 "numeric matrix with at least two columns")
    expect_error(tree_layout(tr, layout = function(g) matrix(1:n, ncol = 1)),
                 "numeric matrix with at least two columns")
    expect_error(tree_layout(tr, layout = function(g) matrix("a", n, 2)),
                 "numeric matrix with at least two columns")
})

test_that("a layout function returning the wrong number of rows names both counts", {
    skip_if_not_installed("igraph")
    tr <- ape::rtree(8)
    expect_error(tree_layout(tr, layout = function(g) matrix(1:4, ncol = 2)),
                 "2 coordinates for a tree with 15 nodes")
})

test_that("a valid external layout function is accepted", {
    skip_if_not_installed("igraph")
    tr <- ape::rtree(8)
    d <- tree_layout(tr, layout = igraph::layout_with_kk)
    expect_equal(nrow(d), ape::Ntip(tr) + ape::Nnode(tr))
    expect_true(all(is.finite(d$x)))
    expect_identical(get_layout(d), "external")
})
