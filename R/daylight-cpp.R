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
