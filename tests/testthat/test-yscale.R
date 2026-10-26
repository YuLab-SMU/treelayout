## `yscale` had no tests. It is a documented parameter with real dispatch --
## numeric column versus category column, with a mapping required for the
## latter -- so a broken branch here would only show up in a plot.

yscale_treedata <- function(tr, data) {
    new("treedata", phylo = tr, data = data)
}

test_that("yscale = 'none' leaves the tip spacing alone", {
    tr <- ape::rtree(8)
    a <- suppressWarnings(tree_layout(tr))
    b <- suppressWarnings(tree_layout(tr, yscale = "none"))
    expect_identical(a$y, b$y)
})

test_that("a numeric yscale column becomes the tip y coordinates", {
    tr <- ape::rtree(8)
    val <- seq(0, 1, length.out = 8)
    td <- yscale_treedata(tr, tibble::tibble(node = 1:8, val = val))
    d <- suppressWarnings(tree_layout(td, yscale = "val"))
    expect_equal(range(d$y[d$isTip]), range(val))
    ## the order should follow the values, not the tip index
    expect_equal(d$label[d$isTip][order(d$y[d$isTip])],
                 tr$tip.label[order(val)])
})

test_that("a category yscale column needs a mapping", {
    tr <- ape::rtree(8)
    td <- yscale_treedata(tr, tibble::tibble(node = 1:8, grp = rep(c("A", "B"), 4)))
    expect_error(suppressWarnings(tree_layout(td, yscale = "grp")),
                 "yscale_mapping")
})

test_that("a category yscale column uses the mapping", {
    tr <- ape::rtree(8)
    td <- yscale_treedata(tr, tibble::tibble(node = 1:8, grp = rep(c("A", "B"), 4)))
    d <- suppressWarnings(tree_layout(td, yscale = "grp",
                                      yscale_mapping = c(A = 1, B = 2)))
    expect_true(all(d$y[d$isTip] %in% c(1, 2)))
})

test_that("a mapping must be named and numeric", {
    tr <- ape::rtree(8)
    td <- yscale_treedata(tr, tibble::tibble(node = 1:8, grp = rep(c("A", "B"), 4)))
    expect_error(suppressWarnings(tree_layout(td, yscale = "grp",
                                              yscale_mapping = 1:2)),
                 "named numeric vector")
})

test_that("an unknown yscale column warns and changes nothing", {
    tr <- ape::rtree(8)
    base <- suppressWarnings(tree_layout(tr))
    expect_warning(d <- tree_layout(tr, yscale = "no_such_column"),
                   "not available in data frame")
    expect_identical(d$y, base$y)
})
