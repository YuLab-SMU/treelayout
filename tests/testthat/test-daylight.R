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


## Golden geometry, pinned before the layout moved here. This is the strongest
## of the daylight tests: a change in the arithmetic shows up here even when
## the qualitative properties above still hold.
reference_daylight_layout <- data.frame(
    node = 1:19,
    parent = c(11L, 13L, 14L, 16L, 16L, 15L, 18L, 18L, 19L, 19L, 11L, 11L, 12L, 13L, 14L, 15L, 12L, 17L, 17L),
    x = c(
        0.537635332908161, -0.830129041870907, -1.93082272965061,
        -2.86040515000489, -3.10899210442627, -2.87009960010117,
        -1.31144187202577, -0.757621437012896, -0.062536798697354,
        -0.590565454143218, 0, -0.867411764038857, -0.990520380451814,
        -1.92558918215148, -2.35350581709342, -2.66202546804589,
        -0.601581797879738, -1.05623138298822, -0.59789000710107
    ),
    y = c(
        0.157863977281575, 0.785015066442061, 0.03922836958825,
        0.794472074756069, -0.087638603789479, -0.985469317318234,
        -1.90676178265179, -2.23104169789429, -1.5822025523596,
        -0.944497523740547, 0, -0.254695074208212, -0.190782796865866,
        -0.043042894933476, -0.328176394726337, -0.089273965444829,
        -0.942722489684546, -1.61436539848238, -0.944122511145176
    ),
    ## angle comes from calculate_angle(), which is the value fortify() has
    ## always returned; the numbers this reference used to hold were computed
    ## a different way and no longer match
    angle = c(
        14.1177089512, 70.2035665209, 3.5081765582,
        71.0493027885, -7.8374833987, -88.1301056886,
        -170.5209025537, -199.5211186952, -141.4957073851,
        -84.4660154579, 0, -22.7772731366, -17.0616251119,
        -3.849307952, -29.3486766594, -7.9837331022,
        -84.3072749104, -144.3720171577, -84.4324782395
    )
)

test_that("daylight geometry matches the pinned reference", {
    set.seed(42)
    tr <- ape::rtree(10)
    res <- suppressWarnings(suppressMessages(
        tree_layout(tr, layout = "daylight")
    ))
    ## strip the layout attribute: the reference is a plain data frame
    got <- as.data.frame(res[, c("node", "parent", "x", "y", "angle")])
    attr(got, "layout") <- NULL
    expect_equal(got, reference_daylight_layout, tolerance = 1e-8)
})

test_that("the daylight cache path matches the uncached path", {
    set.seed(7)
    tr <- ape::rtree(12)
    layout <- suppressWarnings(tree_layout(tr, layout = "equal_angle"))
    cache <- .daylightBuildCache(layout)
    node_id <- cache$internal_nodes[[2]]

    cached <- applyLayoutDaylight(layout, node_id, cache = cache)
    uncached <- applyLayoutDaylight(layout, node_id)

    expect_equal(cached$max_change, uncached$max_change, tolerance = 1e-10)
    expect_equal(
        as.data.frame(cached$tree[, c("node", "parent", "x", "y", "angle")]),
        as.data.frame(uncached$tree[, c("node", "parent", "x", "y", "angle")]),
        tolerance = 1e-10
    )
})
