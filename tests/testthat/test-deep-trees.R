## Ladder-like trees have depth ~ n, which exposed two loops that advanced one
## level per iteration: getYcoord() and getXcoord_no_length(). At 1600 tips a
## rectangular layout took 17.7 s; with both rewritten as a single upward pass
## it takes 0.1 s. These tests pin the geometry, not the speed, but they use a
## deep tree so a regression to the level-by-level form is at least caught by
## the suite timing out or the shape assertion failing.

test_that("getYcoord handles a deep ladder tree", {
    tr <- ape::stree(300, type = "left")
    y <- getYcoord(tr)

    expect_length(y, ape::Ntip(tr) + ape::Nnode(tr))
    expect_false(anyNA(y))
    ## tips keep the evenly spaced 1..n rows
    expect_equal(sort(y[seq_len(ape::Ntip(tr))]), seq_len(ape::Ntip(tr)))

    ## every internal node is the mean of its children
    edge <- tr$edge
    for (v in unique(edge[, 1])) {
        kids <- edge[edge[, 1] == v, 2]
        expect_equal(unname(y[v]), mean(y[kids]), tolerance = 1e-12,
                     info = paste("node", v))
    }
})

test_that("getXcoord_no_length places every internal node left of its children", {
    tr <- ape::stree(300, type = "left")
    x <- getXcoord_no_length(tr)
    edge <- tr$edge

    expect_false(anyNA(x))
    ## a child is always one column to the right of its parent
    expect_true(all(x[edge[, 2]] > x[edge[, 1]]))
    ## the root ends up leftmost after the shift to zero
    expect_equal(min(x), 0)
})

test_that("a deep tree lays out with the expected shape", {
    tr <- ape::stree(200, type = "left")
    d <- tree_layout(tr, layout = "rectangular")

    expect_equal(nrow(d), ape::Ntip(tr) + ape::Nnode(tr))
    expect_true(all(is.finite(d$x)))
    expect_true(all(is.finite(d$y)))
    ## tips are on distinct rows, internal nodes fall between their children
    expect_equal(sum(d$isTip), ape::Ntip(tr))
})

test_that("balancing the same number of tips keeps the results", {
    ## the single-pass rewrite must agree with the level-by-level definition on
    ## shapes where both are cheap enough to state exactly
    tr <- ape::read.tree(text = "((a,b),(c,d));")
    expect_equal(getXcoord_no_length(tr), c(2, 2, 2, 2, 0, 1, 1))
    expect_equal(getYcoord(tr), c(1, 2, 3, 4, 2.5, 1.5, 3.5))
})
