## Compiled version of the inner loop of .daylightGetArcAngles().
##
## Kept as a separate entry point so the R implementation stays the
## reference: the tests compare the two with expect_identical.
.daylight_use_cpp <- function() {
    isTRUE(getOption("treelayout.daylight_cpp", TRUE))
}

## Returns c(left = , right = ), the same shape as the R version.
.daylight_arc_cpp <- function(df, origin_id, subtree, cache,
                              theta_left, theta_right) {
    daylight_arc_angles_cpp(
        df_x = as.numeric(df$x),
        df_y = as.numeric(df$y),
        row_index = as.integer(cache$row_index),
        origin_id = as.integer(origin_id),
        subtree_node_ids = as.integer(subtree$subtree),
        children = cache$children,
        has_children = as.logical(cache$has_children),
        theta_left_start = as.numeric(theta_left),
        theta_right_start = as.numeric(theta_right)
    )
}


## Compiled version of the force-directed relaxation in layoutTreeAndLeaf().
## Uses the same option switch as the daylight fast path.
.treeAndLeafRelax_cpp <- function(positions, cache, params) {
    n2 <- nrow(positions)
    tree_and_leaf_relax_cpp(
        positions = positions,
        depth_norm = as.numeric(cache$depth_norm),
        tip_rows = as.integer(cache$tip_rows),
        internal_rows = as.integer(cache$internal_rows),
        root_row = as.integer(cache$root_row),
        edge_parent_rows = as.integer(cache$edge_parent_rows),
        edge_child_rows = as.integer(cache$edge_child_rows),
        target_edge_length = as.numeric(cache$target_edge_length),
        outward_dir = as.numeric(t(cache$outward_dir)),
        mobility = as.numeric(cache$mobility),
        node_ids = as.integer(cache$node_ids),
        leaf_force = as.numeric(params$leaf_force),
        edge_force = as.numeric(params$edge_force),
        anchor_force = as.numeric(params$anchor_force),
        internal_anchor = as.numeric(params$internal_anchor),
        outward_force = as.numeric(params$outward_force),
        step = as.numeric(params$step),
        cooling = as.numeric(params$cooling),
        tol = as.numeric(params$tol),
        max_iter = as.integer(params$max_iter)
    )
}
