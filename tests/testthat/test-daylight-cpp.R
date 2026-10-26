## The compiled path must give bit-identical results to the R reference.
##
## Off by default in these tests would hide regressions, so both are checked
## explicitly: with the switch on, the output must equal the R output exactly.
## identical(), not all.equal() -- atan2() is a pure function and the loop
## keeps its original order, so nothing is allowed to move.

test_that("compiled daylight matches the R implementation exactly", {
    for (n in c(10, 20, 50, 100)) {
        for (s in 1:3) {
            set.seed(s)
            tr <- ape::rtree(n)

            old <- getOption("treelayout.daylight_cpp")
            on.exit(options(treelayout.daylight_cpp = old), add = TRUE)

            options(treelayout.daylight_cpp = FALSE)
            ref <- suppressWarnings(tree_layout(tr, layout = "daylight"))

            options(treelayout.daylight_cpp = TRUE)
            got <- suppressWarnings(tree_layout(tr, layout = "daylight"))

            expect_identical(got$x, ref$x, info = sprintf("x n=%d seed=%d", n, s))
            expect_identical(got$y, ref$y, info = sprintf("y n=%d seed=%d", n, s))
            expect_identical(got$angle, ref$angle, info = sprintf("angle n=%d seed=%d", n, s))
        }
    }
})

test_that("compiled daylight matches on extreme tree shapes", {
    shapes <- list(
        balanced = ape::stree(32, type = "balanced"),
        left     = ape::stree(32, type = "left"),
        star     = ape::stree(20, type = "star")
    )
    for (nm in names(shapes)) {
        old <- getOption("treelayout.daylight_cpp")
        on.exit(options(treelayout.daylight_cpp = old), add = TRUE)

        options(treelayout.daylight_cpp = FALSE)
        ref <- suppressWarnings(tree_layout(shapes[[nm]], layout = "daylight"))
        options(treelayout.daylight_cpp = TRUE)
        got <- suppressWarnings(tree_layout(shapes[[nm]], layout = "daylight"))

        expect_identical(got$x, ref$x, info = nm)
        expect_identical(got$y, ref$y, info = nm)
    }
})

test_that("compiled daylight respects MAX_COUNT the same way", {
    set.seed(1)
    tr <- ape::rtree(40)
    for (k in c(1, 3, 5)) {
        old <- getOption("treelayout.daylight_cpp")
        on.exit(options(treelayout.daylight_cpp = old), add = TRUE)

        options(treelayout.daylight_cpp = FALSE)
        ref <- suppressWarnings(tree_layout(tr, layout = "daylight", MAX_COUNT = k))
        options(treelayout.daylight_cpp = TRUE)
        got <- suppressWarnings(tree_layout(tr, layout = "daylight", MAX_COUNT = k))

        expect_identical(got$x, ref$x, info = paste("MAX_COUNT", k))
        expect_identical(got$y, ref$y, info = paste("MAX_COUNT", k))
    }
})

test_that("tree_and_leaf inherits the compiled daylight seed layout", {
    ## tree_and_leaf starts from layoutDaylight(); switching the compiled path
    ## must not change its starting point, and therefore not its result.
    set.seed(2)
    tr <- ape::rtree(30)
    old <- getOption("treelayout.daylight_cpp")
    on.exit(options(treelayout.daylight_cpp = old), add = TRUE)

    options(treelayout.daylight_cpp = FALSE)
    ref <- suppressWarnings(tree_layout(tr, layout = "tree_and_leaf"))
    options(treelayout.daylight_cpp = TRUE)
    got <- suppressWarnings(tree_layout(tr, layout = "tree_and_leaf"))

    expect_identical(got$x, ref$x)
    expect_identical(got$y, ref$y)
})
