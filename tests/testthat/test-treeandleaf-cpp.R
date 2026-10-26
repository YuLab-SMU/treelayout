## The compiled relaxation must agree with the R reference within floating
## point noise.
##
## identical() is the wrong bar here, on purpose: the R loop gathers each
## tip's repulsion through colSums() while C++ accumulates pair by pair, so
## addition happens in a different order and the results differ in the last
## bit (~1e-15). 1e-9 is still six orders of magnitude below anything visible,
## and the relaxation has its own convergence tolerance far above it.

test_that("compiled relaxation matches the R implementation", {
    for (n in c(10, 20, 60, 150)) {
        set.seed(1)
        tr <- ape::rtree(n)
        old <- getOption("treelayout.daylight_cpp")
        on.exit(options(treelayout.daylight_cpp = old), add = TRUE)

        options(treelayout.daylight_cpp = FALSE)
        ref <- suppressWarnings(tree_layout(tr, layout = "tree_and_leaf"))
        options(treelayout.daylight_cpp = TRUE)
        got <- suppressWarnings(tree_layout(tr, layout = "tree_and_leaf"))

        expect_equal(got$x, ref$x, tolerance = 1e-9, info = paste("n =", n))
        expect_equal(got$y, ref$y, tolerance = 1e-9, info = paste("n =", n))
    }
})

test_that("compiled relaxation agrees on extreme shapes", {
    shapes <- list(
        balanced = ape::stree(32, type = "balanced"),
        left     = ape::stree(32, type = "left"),
        star     = ape::stree(20, type = "star")
    )
    for (nm in names(shapes)) {
        old <- getOption("treelayout.daylight_cpp")
        on.exit(options(treelayout.daylight_cpp = old), add = TRUE)

        options(treelayout.daylight_cpp = FALSE)
        ref <- suppressWarnings(tree_layout(shapes[[nm]], layout = "tree_and_leaf"))
        options(treelayout.daylight_cpp = TRUE)
        got <- suppressWarnings(tree_layout(shapes[[nm]], layout = "tree_and_leaf"))

        expect_equal(got$x, ref$x, tolerance = 1e-9, info = nm)
        expect_equal(got$y, ref$y, tolerance = 1e-9, info = nm)
    }
})

test_that("compiled relaxation keeps the root pinned", {
    ## the R loop writes initial_positions[root] back every iteration; the C++
    ## port must do the same or the whole layout drifts. The pinned point is
    ## the *seed* layout's root position -- tree_and_leaf starts from daylight,
    ## not equal_angle, so the root is not at the origin.
    set.seed(3)
    tr <- ape::rtree(40)
    old <- getOption("treelayout.daylight_cpp")
    on.exit(options(treelayout.daylight_cpp = old), add = TRUE)

    seed <- suppressWarnings(tree_layout(tr, layout = "daylight"))
    root <- seed$node[seed$node == seed$parent][1]

    options(treelayout.daylight_cpp = TRUE)
    d <- suppressWarnings(tree_layout(tr, layout = "tree_and_leaf"))

    expect_equal(d$x[d$node == root], seed$x[seed$node == root], tolerance = 1e-12)
    expect_equal(d$y[d$node == root], seed$y[seed$node == root], tolerance = 1e-12)
})

test_that("compiled relaxation is deterministic", {
    set.seed(4)
    tr <- ape::rtree(50)
    old <- getOption("treelayout.daylight_cpp")
    on.exit(options(treelayout.daylight_cpp = old), add = TRUE)

    options(treelayout.daylight_cpp = TRUE)
    a <- suppressWarnings(tree_layout(tr, layout = "tree_and_leaf"))
    b <- suppressWarnings(tree_layout(tr, layout = "tree_and_leaf"))
    expect_identical(a$x, b$x)
    expect_identical(a$y, b$y)
})
