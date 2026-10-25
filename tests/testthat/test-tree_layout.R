test_that("every layout name produces a full node table", {
    tr <- ape::rtree(8, tip.label = paste0("t", 1:8))
    n_expected <- 8 + ape::Nnode(tr)

    for (ly in treelayout:::.layout_names$all) {
        res <- suppressWarnings(tree_layout(tr, layout = ly))
        expect_equal(nrow(res), n_expected, info = ly)
        expect_true(all(c("node", "parent", "x", "y") %in% names(res)), info = ly)
        expect_identical(attr(res, "layout"), ly)
    }
})

test_that("coordinates are finite for every layout", {
    tr <- ape::rtree(10)
    for (ly in treelayout:::.layout_names$all) {
        res <- suppressWarnings(tree_layout(tr, layout = ly))
        expect_true(all(is.finite(res$x)), info = ly)
        expect_true(all(is.finite(res$y)), info = ly)
    }
})

test_that("an invalid layout name is rejected", {
    ## before the package existed, fortify() returned silently for this
    expect_error(tree_layout(ape::rtree(5), layout = "totally_fake_layout"))
})

test_that("branch.length = 'none' gives a cladogram", {
    tr <- ape::rtree(6)
    a <- tree_layout(tr, branch.length = "none")
    b <- tree_layout(tr, branch.length = "branch.length")
    expect_false(isTRUE(all.equal(a$x, b$x)))
})

test_that("ladderize toggles tip order", {
    tr <- ape::rtree(8)
    a <- tree_layout(tr, ladderize = TRUE)
    b <- tree_layout(tr, ladderize = FALSE)
    expect_equal(nrow(a), nrow(b))
})

test_that("treedata annotations are carried through", {
    ## as.treedata.phylo lives in treeio, so build the object directly
    tr <- ape::rtree(5)
    td <- tidytree::treedata(phylo = tr,
                             data = tibble::tibble(
                                 node = seq_len(5 + ape::Nnode(tr)),
                                 value = seq_len(5 + ape::Nnode(tr))))
    res <- tree_layout(td)
    expect_true("value" %in% names(res))
})

test_that("get_layout reads the layout name", {
    tr <- ape::rtree(6)
    res <- tree_layout(tr, layout = "circular")
    expect_identical(get_layout(res), "circular")
    ## also works on anything with a $data slot
    fake <- list(data = res)
    expect_identical(get_layout(fake), "circular")
})

test_that("unrooted grouping matches the old hard-coded list", {
    expect_setequal(.layout_names$unrooted,
                    c("equal_angle", "daylight", "ape", "tree_and_leaf"))
})

test_that("package has no plotting or Bioconductor dependency", {
    ## read DESCRIPTION directly: `desc` would add a dependency just for a test
    d <- read.dcf(system.file("DESCRIPTION", package = "treelayout"))
    imports <- strsplit(paste(d[1, c("Imports", "Depends")], collapse = ","), ",")[[1]]
    imports <- trimws(sub("\\(.*", "", imports))
    banned <- c("ggplot2", "grid", "aplot", "ggfun", "ggtree", "ggstar",
                "treeio", "Biostrings", "BiocGenerics", "S4Vectors")
    expect_false(any(banned %in% imports))
})

test_that("multiPhylo is laid out tree by tree and stacked", {
    trs <- ape::rmtree(3, 5)
    res <- tree_layout(trs)
    n_expected <- 3 * (5 + ape::Nnode(trs[[1]]))
    expect_equal(nrow(res), n_expected)
    expect_true(".id" %in% names(res))
    expect_equal(nlevels(res$.id), 3)
    expect_identical(attr(res, "layout"), "rectangular")
})

test_that("named multiPhylo keeps the names", {
    trs <- ape::rmtree(2, 4)
    names(trs) <- c("a", "b")
    res <- tree_layout(trs)
    expect_identical(levels(res$.id), c("a", "b"))
})

test_that("data can be omitted", {
    expect_silent(tree_layout(ape::rtree(5)))
})

test_that("daylight layout is quiet unless the verbose option is set", {
    tr <- ape::rtree(8)
    msgs <- character()
    withCallingHandlers(
        invisible(tree_layout(tr, layout = "daylight")),
        message = function(m) {
            msgs <<- c(msgs, conditionMessage(m))
            invokeRestart("muffleMessage")
        })
    expect_length(msgs, 0)
})

test_that("multiPhylo result has a clean class vector", {
    res <- tree_layout(ape::rmtree(3, 4))
    expect_identical(class(res), c("tbl_tree", "tbl_df", "tbl", "data.frame"))
})
