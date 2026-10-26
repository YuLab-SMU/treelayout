## Frozen copy of .daylightGetArcAngles() before vectorisation.
##
## The vectorised version must return bit-identical angles, and that can only
## be checked against the original. This file is loaded only by the tests
## (tests/testthat/helper-*.R is sourced automatically) and is not part of
## the package.
##
## Do not edit: it is a reference, not an implementation.

.daylightGetArcAngles_reference <- function(df, origin_id, subtree, cache) {
  df_x <- df$x
  df_y <- df$y
  row_index <- cache$row_index
  origin_row <- row_index[origin_id]
  x_origin <- df_x[origin_row]
  y_origin <- df_y[origin_row]
  subtree_root_id <- subtree$node
  subtree_node_ids <- subtree$subtree
  children_ids <- cache$children[[origin_id]]

  theta_left <- NA_real_
  theta_right <- NA_real_

  if (subtree_root_id %in% children_ids) {
    subtree_root_row <- row_index[subtree_root_id]
    theta_left <- getNodeAngle.vector(x_origin, y_origin, df_x[subtree_root_row], df_y[subtree_root_row])
    theta_right <- theta_left
  } else if (subtree_root_id == origin_id) {
    if (length(children_ids) == 0L) {
      return(c(left = 0, right = 0))
    }
    if (length(children_ids) == 2L) {
      child_rows <- row_index[children_ids]
      theta1 <- getNodeAngle.vector(x_origin, y_origin, df_x[child_rows[1]], df_y[child_rows[1]])
      theta2 <- getNodeAngle.vector(x_origin, y_origin, df_x[child_rows[2]], df_y[child_rows[2]])
      delta <- theta1 - theta2
      if (delta > 1) {
        delta_adj <- delta - 2
      } else if (delta < -1) {
        delta_adj <- delta + 2
      } else {
        delta_adj <- delta
      }
      if (delta_adj >= 0) {
        theta_left <- theta1
        theta_right <- theta2
      } else {
        theta_left <- theta2
        theta_right <- theta1
      }
    } else {
      child_row <- row_index[children_ids[1]]
      theta_left <- getNodeAngle.vector(x_origin, y_origin, df_x[child_row], df_y[child_row])
      theta_right <- theta_left
    }
  } else {
    tree_root <- cache$root_node
    if (!is.na(tree_root)) {
      tree_root_row <- row_index[tree_root]
      theta_left <- getNodeAngle.vector(x_origin, y_origin, df_x[tree_root_row], df_y[tree_root_row])
      theta_right <- theta_left
    } else {
      return(c(left = 0, right = 0))
    }
  }

  if (is.na(theta_left) || length(subtree_node_ids) == 0L) {
    return(c(left = 0, right = 0))
  }

  arc <- c(left = theta_left, right = theta_right)
  parent_nodes <- subtree_node_ids[subtree_node_ids != origin_id]
  parent_nodes <- parent_nodes[cache$has_children[parent_nodes]]

  for (parent_id in parent_nodes) {
    parent_row <- row_index[parent_id]
    theta_parent <- getNodeAngle.vector(x_origin, y_origin, df_x[parent_row], df_y[parent_row])
    child_ids <- cache$children[[parent_id]]
    child_ids <- child_ids[child_ids != origin_id]
    if (length(child_ids) == 0L) {
      next
    }

    for (child_id in child_ids) {
      child_row <- row_index[child_id]
      theta_child <- getNodeAngle.vector(x_origin, y_origin, df_x[child_row], df_y[child_row])
      if ((arc['left'] < arc['right'] && !(theta_child > arc['left'] && theta_child < arc['right'])) ||
          (arc['left'] > arc['right'] && (theta_child < arc['left'] && theta_child > arc['right']))) {
        next
      }

      delta <- theta_child - theta_parent
      delta_adj <- delta
      if (delta > 1) {
        delta_adj <- delta - 2
      } else if (delta < -1) {
        delta_adj <- delta + 2
      }

      theta_child_adj <- theta_child
      if (delta_adj > 0) {
        if (abs(delta) > 1) {
          if (arc['left'] > 0 && theta_child < 0) {
            theta_child_adj <- theta_child + 2
          } else if (arc['left'] < 0 && theta_child > 0) {
            theta_child_adj <- theta_child - 2
          }
        }
        if (arc['left'] < theta_child_adj) {
          arc['left'] <- theta_child
        }
      } else if (delta_adj < 0) {
        if (abs(delta) > 1) {
          if (arc['right'] > 0 && theta_child < 0) {
            theta_child_adj <- theta_child + 2
          } else if (arc['right'] < 0 && theta_child > 0) {
            theta_child_adj <- theta_child - 2
          }
        }
        if (arc['right'] > theta_child_adj) {
          arc['right'] <- theta_child
        }
      }
    }
  }

  arc[arc < 0] <- arc[arc < 0] + 2
  arc
}
