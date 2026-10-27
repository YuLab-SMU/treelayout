## ggtree#400 was "C stack usage is too close to the limit" with layout =
## "daylight". The cause was a mutate.tbl_tree() that dispatched back into
## itself in an old tidytree, so it was never about the layout code. A later
## comment on the same issue reported it again on a 3400-tip tree with
## geom_inset, which is the plotting side.
##
## These tests do not prove the absence of recursion, but they do fail loudly
## if a layout starts recursing per tree level: a ladder tree has depth ~ n, so
## anything walking it recursively blows the C stack long before n = 3000.

test_that("a deep ladder tree does not exhaust the C stack", {
    tr <- ape::stree(3000, type = "left")
    for (ly in c("rectangular", "circular", "fan", "equal_angle", "slanted", "tidy")) {
        res <- tryCatch(suppressWarnings(tree_layout(tr, layout = ly)),
                        error = function(e) conditionMessage(e))
        expect_false(is.character(res), info = ly)
        expect_equal(nrow(res), ape::Ntip(tr) + ape::Nnode(tr), info = ly)
    }
})

test_that("the tidytree helpers the layouts rely on are iterative", {
    ## offspring_bulk() is the one that used to recurse per level; it now walks
    ## a post-order with a for loop. Guard the property we depend on.
    tr <- ape::stree(2000, type = "left")
    expect_no_error(offspring(tr, ape::Ntip(tr) + 1L, self_include = TRUE))
})

test_that("a wide tree with many tips does not exhaust the C stack", {
    tr <- ape::rtree(3400)
    expect_no_error(suppressWarnings(tree_layout(tr, layout = "circular")))
})
