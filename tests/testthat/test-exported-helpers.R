## The exported helpers had no tests of their own -- they were only reached
## indirectly through tree_layout(), or not at all. Each one is public API, so
## it needs a contract test that fails if its behaviour changes.

test_that("Date2decimal and decimal2Date round-trip", {
    d <- as.Date(c("2014-05-05", "2020-01-01", "1999-12-31"))
    dec <- Date2decimal(d)
    expect_type(dec, "double")
    expect_equal(dec, c(2014.34, 2020, 1999.997), tolerance = 1e-3)
    expect_equal(decimal2Date(dec), d)
})

test_that("Date2decimal passes numeric through unchanged", {
    ## documented: a numeric is returned as is, so callers can apply it to
    ## values that are already decimal
    expect_identical(Date2decimal(2014.34), 2014.34)
})

test_that("decimal2Date handles the first day of a year", {
    expect_equal(decimal2Date(2014), as.Date("2014-01-01"))
})

test_that("calculate_angle adds an angle column", {
    d <- suppressWarnings(tree_layout(ape::rtree(10), layout = "circular"))
    out <- calculate_angle(d)
    expect_true("angle" %in% names(out))
    expect_equal(nrow(out), nrow(d))
    expect_true(all(is.finite(out$angle)))
})

test_that("calculate_angle rejects data without a y column", {
    d <- suppressWarnings(tree_layout(ape::rtree(6), layout = "circular"))
    expect_error(calculate_angle(d[, setdiff(names(d), "y")]), "y")
})

test_that("calculate_branch_mid adds branch only for aligned layouts", {
    tr <- ape::rtree(10)
    for (ly in c("circular", "rectangular", "dendrogram")) {
        d <- suppressWarnings(tree_layout(tr, layout = ly))
        out <- calculate_branch_mid(d, ly)
        expect_true("branch" %in% names(out), info = ly)
        expect_false("branch.x" %in% names(out), info = ly)
        expect_false("branch.y" %in% names(out), info = ly)
    }
})

test_that("calculate_branch_mid adds branch.x and branch.y for unrooted layouts", {
    ## the layout name list drives this; if someone edits it, the columns
    ## silently appear or vanish, so pin the behaviour
    tr <- ape::rtree(10)
    for (ly in c("equal_angle", "daylight", "ape", "tree_and_leaf")) {
        d <- suppressWarnings(tree_layout(tr, layout = ly))
        out <- calculate_branch_mid(d, ly)
        expect_true("branch" %in% names(out), info = ly)
        expect_true("branch.x" %in% names(out), info = ly)
        expect_true("branch.y" %in% names(out), info = ly)
    }
})

test_that("re_assign_ycoord_df fills missing y bottom-up", {
    d <- suppressWarnings(tree_layout(ape::rtree(10), layout = "rectangular"))
    ## blank out an internal node's y and let the helper put it back
    internal <- which(!d$isTip)[1]
    expected <- d$y[internal]
    d$y[internal] <- NA_real_
    out <- re_assign_ycoord_df(d, currentNode = d$node[d$isTip])
    expect_false(anyNA(out$y))
    expect_equal(out$y[internal], expected, tolerance = 1e-12)
})

test_that("re_assign_ycoord_df leaves present y values alone", {
    d <- suppressWarnings(tree_layout(ape::rtree(8), layout = "rectangular"))
    out <- re_assign_ycoord_df(d, currentNode = d$node[d$isTip])
    expect_equal(out$y, d$y)
})

test_that("getYcoord places tips on evenly spaced rows", {
    tr <- ape::rtree(10)
    y <- getYcoord(tr)
    expect_length(y, ape::Ntip(tr) + ape::Nnode(tr))
    expect_equal(sort(y[seq_len(ape::Ntip(tr))]), seq_len(ape::Ntip(tr)))
})

test_that("getYcoord honours tip.order", {
    tr <- ape::rtree(10)
    ord <- rev(tr$tip.label)
    y <- getYcoord(tr, tip.order = ord)
    ## the first tip label in `ord` should get the smallest y
    expect_equal(unname(y[match(ord[1], tr$tip.label)]), 1)
})

test_that("getYcoord fills internal nodes with the mean of their children", {
    tr <- ape::rtree(8)
    y <- getYcoord(tr)
    edge <- tr$edge
    for (v in unique(edge[, 1])) {
        kids <- edge[edge[, 1] == v, 2]
        expect_equal(unname(y[v]), mean(y[kids]), tolerance = 1e-12,
                     info = paste("node", v))
    }
})

test_that("getXcoord2 walks the edge list from the root", {
    tr <- ape::rtree(8)
    x <- numeric(max(tr$edge))
    out <- getXcoord2(x, root = ape::Ntip(tr) + 1L,
                      parent = tr$edge[, 1], child = tr$edge[, 2],
                      len = tr$edge.length)
    expect_length(out, max(tr$edge))
    expect_false(anyNA(out))
    expect_equal(out[ape::Ntip(tr) + 1L], 0)
})

test_that("getXcoord2 respects the start argument", {
    tr <- ape::rtree(6)
    x <- numeric(max(tr$edge))
    out <- getXcoord2(x, root = ape::Ntip(tr) + 1L,
                      parent = tr$edge[, 1], child = tr$edge[, 2],
                      len = tr$edge.length, start = 5)
    expect_equal(out[ape::Ntip(tr) + 1L], 5)
})
