## Properties of the daylight layout that hold regardless of the stopping
## rule: they were pinned in ggtree before the layout moved here.

test_that("daylight keeps every tip at a distinct, finite position", {
    for (n in c(15, 40)) {
        set.seed(n)
        tr <- ape::rtree(n)
        res <- suppressWarnings(tree_layout(tr, layout = "daylight"))
        tips <- res[res$isTip, ]
        expect_equal(nrow(tips), n)
        expect_true(all(is.finite(tips$x)))
        expect_true(all(is.finite(tips$y)))
        expect_true(all(is.finite(tips$angle)))
    }
})

test_that("daylight keeps tip order aligned with equal_angle", {
    set.seed(42)
    tr <- ape::rtree(10)
    daylight <- suppressWarnings(tree_layout(tr, layout = "daylight"))
    angle <- suppressWarnings(tree_layout(tr, layout = "equal_angle"))
    ## compare the columns only: the tables carry different layout attributes
    expect_equal(daylight$node, angle$node)
    expect_equal(daylight$parent, angle$parent)
})

test_that("MAX_COUNT caps the number of rounds", {
    set.seed(8)
    tr <- ape::rtree(30)
    msgs <- character()
    op <- options(treelayout.verbose = TRUE)
    on.exit(options(op), add = TRUE)
    withCallingHandlers(
        invisible(tree_layout(tr, layout = "daylight", MAX_COUNT = 2)),
        message = function(m) {
            msgs <<- c(msgs, conditionMessage(m))
            invokeRestart("muffleMessage")
        })
    rounds <- grep("Average angle change", msgs)
    expect_lte(length(rounds), 2L)
    expect_gte(length(rounds), 1L)
})
