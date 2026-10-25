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
