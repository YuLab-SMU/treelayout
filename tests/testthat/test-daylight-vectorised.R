## The vectorised .daylightGetArcAngles() must return the same angles as the
## original. helper-daylight-reference.R holds a frozen copy of the original
## for that comparison.
##
## Bit-identical, not approximate: atan2() is a pure function, so vectorising
## it must not move a single bit. If this ever needs all.equal, the change is
## not a pure vectorisation and has to be justified separately.

## Build the cache and a few subtrees the way layoutDaylight() does, so the
## test exercises the real call shape rather than a hand-made one.
daylight_case <- function(n, seed) {
    set.seed(seed)
    df <- suppressWarnings(tree_layout(ape::rtree(n), layout = "equal_angle"))
    cache <- .daylightBuildCache(df)
    list(df = df, cache = cache, internal = cache$internal_nodes)
}

test_that("vectorised arc angles match the reference on random trees", {
    for (n in c(10, 20, 50, 100)) {
        for (s in 1:3) {
            cs <- daylight_case(n, s)
            checked <- 0L
            for (origin in cs$internal) {
                subs <- cs$cache$subtrees[[origin]]
                if (is.null(subs)) next
                for (st in subs) {
                    old <- .daylightGetArcAngles_reference(cs$df, origin, st, cs$cache)
                    new <- .daylightGetArcAngles(cs$df, origin, st, cs$cache)
                    expect_identical(new, old,
                                     info = sprintf("n=%d seed=%d origin=%d", n, s, origin))
                    checked <- checked + 1L
                }
            }
            expect_gt(checked, 0L)
        }
    }
})

test_that("vectorised arc angles match the reference on extreme tree shapes", {
    shapes <- list(
        balanced = ape::stree(32, type = "balanced"),
        left     = ape::stree(32, type = "left"),
        star     = ape::stree(20, type = "star")
    )
    for (nm in names(shapes)) {
        df <- suppressWarnings(tree_layout(shapes[[nm]], layout = "equal_angle"))
        cache <- .daylightBuildCache(df)
        for (origin in cache$internal_nodes) {
            subs <- cache$subtrees[[origin]]
            if (is.null(subs)) next
            for (st in subs) {
                expect_identical(
                    .daylightGetArcAngles(df, origin, st, cache),
                    .daylightGetArcAngles_reference(df, origin, st, cache),
                    info = paste(nm, origin)
                )
            }
        }
    }
})

test_that("vectorised arc angles match on the sizes where the old early stop drifted", {
    ## the reverted early stop changed results by ~99% of the x range at
    ## n = 200 and ~101% at n = 300, so those sizes go in the test set
    for (n in c(200, 300)) {
        cs <- daylight_case(n, 1)
        for (origin in cs$internal) {
            subs <- cs$cache$subtrees[[origin]]
            if (is.null(subs)) next
            for (st in subs) {
                expect_identical(
                    .daylightGetArcAngles(cs$df, origin, st, cs$cache),
                    .daylightGetArcAngles_reference(cs$df, origin, st, cs$cache),
                    info = sprintf("n=%d origin=%d", n, origin)
                )
            }
        }
    }
})
