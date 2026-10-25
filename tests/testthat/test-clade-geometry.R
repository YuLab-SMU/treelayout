## helper: the root node of a node table
treelayout_root <- function(data) {
    data$node[data$node == data$parent][1]
}

## clade_geometry() replaces three exported functions
## (getTreeArcAngles / getNodeEuclDistances / .getRoot.df). The comparison
## below calls those three directly: they run in this package's namespace,
## so the internals are reachable even after they stop being exported.

test_that("clade_geometry() matches the three functions it replaces", {
    set.seed(1)
    tr <- ape::rtree(15, tip.label = paste0("t", 1:15))
    data <- suppressWarnings(tree_layout(tr, layout = "circular"))

    ## a couple of internal nodes and the root
    internal <- unique(data$parent[data$parent != data$node])
    nodes <- c(internal[1], internal[2], treelayout_root(data))

    for (node in nodes) {
        got <- clade_geometry(data, node)

        subtree <- list(subtree = getSubtree.df(data, node), node = node)
        arc <- getTreeArcAngles(data, node, subtree)

        expect_equal(got$left, as.numeric(arc["left"]), info = paste("left", node))
        expect_equal(got$right, as.numeric(arc["right"]), info = paste("right", node))
        expect_equal(
            got$radius,
            max(getNodeEuclDistances(data[data$node %in% subtree$subtree, ], node)),
            info = paste("radius", node)
        )
        expect_equal(
            got$root_radius,
            max(getNodeEuclDistances(data, .getRoot.df(data))),
            info = paste("root_radius", node)
        )
    }
})
