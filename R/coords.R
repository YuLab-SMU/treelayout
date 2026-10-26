

##' @importFrom ape reorder.phylo
layout.unrooted <- function(model, branch.length="branch.length", layout.method="equal_angle", MAX_COUNT=5, ...) {

    df <- switch(layout.method,
                 equal_angle = layoutEqualAngle(model, branch.length),
                 daylight = layoutDaylight(model, branch.length, MAX_COUNT),
                 tree_and_leaf = layoutTreeAndLeaf(model, branch.length, MAX_COUNT = MAX_COUNT, ...),
    			 ape = layoutApe(model, branch.length)
                 )

    return(df)
}

set_branch_length_cladogram <- function(tree) {
    phylo <- as.phylo(tree)
    edge <- phylo$edge
    xpos <- getXcoord_no_length(phylo)
    phylo$edge.length <- xpos[edge[,2]] - xpos[edge[,1]]

    if (is(tree, "phylo")) {
        return(phylo)
    } else if (is(tree, "treedata")) {
        tree@phylo <- phylo
        return(tree)
    }
    message("unknown tree object, fail to set branch length for cladogram...")
    return(tree)
}

## Node -> children and node -> parent as plain lookup structures, built once.
##
## `child(tree, v)` / `parent(tree, v)` scan the whole edge list on every call,
## so using them inside a per-node loop is O(n^2) on a ladder-like tree. These
## give the same answers in O(1) per lookup.
.children_index <- function(tree) {
    edge <- tree[["edge"]]
    parent <- edge[, 1]
    child <- edge[, 2]
    n <- max(edge)
    cl <- split(child, parent)
    out <- vector("list", n)
    out[as.integer(names(cl))] <- cl
    for (i in seq_len(n)) {
        out[[i]] <- if (is.null(out[[i]])) integer() else as.integer(out[[i]])
    }
    out
}

.parent_index <- function(tree) {
    edge <- tree[["edge"]]
    n <- max(edge)
    out <- integer(n)
    out[edge[, 2]] <- edge[, 1]
    out
}

## Slanted coordinates straight from the tree.
##
## The generic path built an n_tip x n_node indicator matrix (sbp) and then took
## colSums() and a weighted colMeans() of it. That matrix is O(n^2) in memory and
## the rbind() to assemble it dominated the profile (74% in `c()`), which is why
## `slanted` took 4.7 s at 2000 ladder tips. Both summaries are simple enough to
## accumulate in one upward pass:
##
##   x[node] = -(number of tips under node)
##   y[node] = mean of the 1-based tip rows under node
##
## Row order matches sbp, whose rows are the tips in `.nodeId(tree, "tips")`
## order; for a phylo that is 1..n_tip, and the tip label order of the tree is
## what getXcoord_no_length_slanted() then re-sorts by node id anyway.
.slanted_coords <- function(tree) {
    edge <- tree[["edge"]]
    parent <- edge[, 1]
    child <- edge[, 2]
    n <- max(edge)
    n_tip <- length(tree[["tip.label"]])

    n_desc <- integer(n)
    sum_row <- numeric(n)
    ## Tips are numbered by their position in the edge list. That is what
    ## `.convert_tips2ancestors_sbp()` did too -- its rows are the tips in
    ## `.nodeId(tree, "tips")` order, which is this same order -- so the two
    ## agree even after reorder(), when the tips no longer run 1..n.
    tip_ids <- child[!(child %in% parent)]
    n_desc[tip_ids] <- 1L
    sum_row[tip_ids] <- seq_along(tip_ids)

    parent_of <- .parent_index(tree)
    for (v in order(ape::node.depth(tree))) {
        p <- parent_of[v]
        if (p != 0L) {
            n_desc[p] <- n_desc[p] + n_desc[v]
            sum_row[p] <- sum_row[p] + sum_row[v]
        }
    }

    xs <- -n_desc
    xs <- xs + max(abs(xs))

    ys <- ifelse(n_desc > 0L, sum_row / pmax(n_desc, 1L), NA_real_)
    ys[n_desc == 0L] <- NA_real_
    ys[is.nan(ys)] <- NA_real_

    list(x = xs, y = ys)
}

## Number of tips in the subtree of every node, indexed by node id.
##
## `layoutEqualAngle()` computes this with
##   vapply(1:N, function(i) length(offspring(tree, i, type = "tips",
##                                           self_include = TRUE)), numeric(1))
## which walks each subtree from scratch and is O(n^2) in the number of nodes:
## 0.34 s at 500 tips, 4.8 s at 2000 tips.  One upward pass over the tree in
## order of increasing depth gives the same counts in O(n).
##
## Two details that are easy to get wrong here:
##  - node ids are not assumed to be `1:N` in tree order, so tips are the
##    child ids that never appear as a parent;
##  - `offspring(..., self_include = TRUE)` counts the node itself, so every
##    internal node gets an extra +1 on top of its descendants.
.count_tips_postorder <- function(tree) {
    edge <- tree$edge
    parent <- edge[, 1]
    child <- edge[, 2]
    n <- max(edge)

    is_tip_edge <- !(child %in% parent)
    tip_ids <- child[is_tip_edge]

    nb <- integer(n)
    nb[tip_ids] <- 1L

    ## node -> its parent, as a lookup vector (0 for the root)
    pvec <- integer(n)
    pvec[child] <- parent

    ## `ape::node.depth()` grows with distance from the root, so increasing
    ## depth visits the deepest nodes first: every child's count is final
    ## before it is added into its parent.
    for (v in order(ape::node.depth(tree))) {
        p <- pvec[v]
        if (p != 0L) {
            nb[p] <- nb[p] + nb[v]
        }
    }

    ## count each internal node itself, matching `self_include = TRUE`
    internal <- which(!(seq_len(n) %in% tip_ids))
    nb[internal] <- nb[internal] + 1L

    nb
}

##' 'Equal-angle layout algorithm for unrooted trees'
##'
##' @references
##' "Inferring Phylogenies" by Joseph Felsenstein.
##'
##' @title layoutEqualAngle
##' @param model tree object, e.g. phylo or treedata
##' @param branch.length set to 'none' for edge length of 1. Otherwise the phylogenetic tree edge length is used.
##' @return tree as data.frame with equal angle layout.
##' @noRd
layoutEqualAngle <- function(model, branch.length = "branch.length"){
	tree <- as.phylo(model)

  if (! is.null(tree$edge.length)) {
      if (anyNA(tree$edge.length)) {
          warning("'edge.length' contains NA values...\n",
                  "## setting 'edge.length' of the tree to NULL ",
                  "automatically when plotting the tree...")
          tree$edge.length <- NULL
      }
  }

  if (is.null(tree$edge.length) || branch.length == "none") {
      tree <- set_branch_length_cladogram(tree)
  }
  N <- getNodeNum(tree)
  brlen <- numeric(N)
  brlen[tree$edge[,2]] <- tree$edge.length

  root <- rootnode(tree)
  ## Convert Phylo tree to data.frame.

  df <- as_tibble(model) |>
      mutate(isTip = ! .data$node %in% .data$parent)

    ## NOTE: Angles (start, end, angle) are in half-rotation units (radians/pi or degrees/180)

    ## The geometry is accumulated in plain numeric vectors rather than in
    ## `df`: assigning into a tbl_df row by row costs a method dispatch
    ## (`[<-.tbl_df` -> `tbl_subassign`) and a data-frame subset for every
    ## single value, which dominated this function (65% of its time in
    ## `Ops.data.frame`).  The vectors are written back once at the end.
    x <- numeric(N)
    y <- numeric(N)
    start <- numeric(N)
    end   <- numeric(N)
    angle <- numeric(N)
    ## Initialize root node position and angles.
    end[root] <- 2 # 360-degrees

    brlen <- brlen[df$node] # for cladogram

    ## Get number of tips for each node in tree.
    nb.sp <- .count_tips_postorder(tree)
    ## Get list of node id's.
    nodes <- getNodes_by_postorder(tree)

    ## children of every node, so the inner loop does not rebuild them
    kids <- child(tree, nodes)

    for (k in seq_along(nodes)) {
        curNode <- nodes[k]
        ## Get number of tips for current node.
        curNtip <- nb.sp[curNode]
        ## Get array of child node indexes of current node.
        children <- kids[[k]]
        if (length(children) == 0L) next

        ## Get "start" and "end" angles of a segment for current node.
        seg_start <- start[curNode]
        total_angle <- end[curNode] - seg_start
        cur_x <- x[curNode]
        cur_y <- y[curNode]

        ## vectorised over the children of this node
        alpha <- total_angle * nb.sp[children] / curNtip
        ## cumulative start of each child segment within the parent's span
        seg_offset <- c(0, cumsum(alpha))[seq_along(alpha)]
        ## beta = angle of line from parent node to i-th child
        beta <- seg_start + seg_offset + alpha / 2
        length.child <- brlen[children]

        x[children] <- cur_x + cospi(beta) * length.child
        y[children] <- cur_y + sinpi(beta) * length.child
        ## orthogonal angle to beta for tip label
        angle[children] <- -90 - 180 * beta * sign(beta - 1)
        ## start and end angles of the child's segment
        start[children] <- seg_start + seg_offset
        end[children]   <- seg_start + cumsum(alpha)
    }

    df$x <- x
    df$y <- y
    df$start <- start
    df$end <- end
    df$angle <- angle
    df$branch.length <- brlen

  tree_df <- as_tibble(df)
  class(tree_df) <- c("tbl_tree", class(tree_df))
  return(tree_df)
}

##' Equal daylight layout method for unrooted trees.
##'
##' #' @title
##' @param model tree object, e.g. phylo or treedata
##' @param branch.length set to 'none' for edge length of 1. Otherwise the phylogenetic tree edge length is used.
##' @param MAX_COUNT the maximum number of iterations to run (default 5)
##' @return tree as data.frame with equal angle layout.
##' @references
##' The following aglorithm aims to implement the vague description of the "Equal-daylight Algorithm"
##' in "Inferring Phylogenies" pp 582-584 by Joseph Felsenstein.
##'
##' ```
##' Leafs are subtrees with no children
##' Initialise tree using equal angle algorithm
##' tree_df = equal_angle(tree)
##'
##' nodes = get list of nodes in tree_df breadth-first
##' nodes = remove tip nodes.
##'
##' ```
##' @noRd
layoutDaylight <- function(model, branch.length, MAX_COUNT=5 ){
	tree <- as.phylo(model)

    ## How to set optimal
    MINIMUM_AVERAGE_ANGLE_CHANGE <- 0.05

    ## Initialize tree.
    tree_df <- layoutEqualAngle(model, branch.length)
    ## The rotation below assigns into df$x / df$y thousands of times. On a
    ## tbl_df each of those goes through tbl_subassign(), which showed up as
    ## ~25% of the profile; a plain data.frame is the same object underneath and
    ## is converted back once the rounds finish.
    ## Keep the tbl_tree class: the cache builder calls offspring(), which
    ## dispatches on it. Only the tbl_df wrapping is dropped, which is what took
    ## the tbl_subassign() detour on every assignment.
    tree_df <- as.data.frame(tree_df, stringsAsFactors = FALSE)
    class(tree_df) <- c("tbl_tree", "data.frame")
    cache <- .daylightBuildCache(tree_df)
    nodes <- cache$internal_nodes

    if (length(nodes) == 0L) {
        tree_df <- as_tibble(tree_df)
        class(tree_df) <- c("tbl_tree", class(tree_df))
        return(tree_df)
    }

    ## Run up to MAX_COUNT rounds, stopping when the average angle change
    ## drops below MINIMUM_AVERAGE_ANGLE_CHANGE.
    ave_change <- 1.0
    for (i in seq_len(MAX_COUNT)) {
        total_max <- 0.0
        for (currentNode_id in nodes) {
            result <- applyLayoutDaylight(tree_df, currentNode_id, cache = cache)
            tree_df <- result$tree
            total_max <- total_max + result$max_change
        }
        ave_change <- total_max / length(nodes)
        ## leftover progress output from the original implementation; keep it
        ## available for debugging but quiet by default
        if (isTRUE(getOption("treelayout.verbose"))) {
            message("Average angle change [", i, "] ", ave_change)
        }
        if (ave_change <= MINIMUM_AVERAGE_ANGLE_CHANGE) break
    }

    tree_df <- as_tibble(tree_df)
    class(tree_df) <- c("tbl_tree", class(tree_df))
    return(tree_df)
}

##' TreeAndLeaf-inspired unrooted layout.
##'
##' @param model tree object, e.g. phylo or treedata
##' @param branch.length set to 'none' for edge length of 1. Otherwise the phylogenetic tree edge length is used.
##' @param MAX_COUNT maximum number of iterations used by the daylight initializer.
##' @param initial_layout initial unrooted layout used before leaf-centric relaxation.
##' @param max_iter maximum number of relaxation iterations.
##' @param leaf_force pairwise repulsion strength among leaves.
##' @param edge_force spring strength used to preserve the tree skeleton.
##' @param anchor_force anchor strength pulling tips toward the initial layout.
##' @param internal_anchor anchor strength pulling internal nodes toward the initial layout.
##' @param outward_force outward radial bias applied to leaves.
##' @param step initial step size.
##' @param cooling multiplicative decay applied to step size per iteration.
##' @param tol convergence tolerance on maximum displacement.
##' @param ... Additional arguments reserved for compatibility with other layout helpers.
##' @return tree as data.frame with a TreeAndLeaf-inspired layout.
##' @noRd
layoutTreeAndLeaf <- function(model,
                              branch.length,
                              MAX_COUNT = 5,
                              initial_layout = "daylight",
                              max_iter = 200L,
                              leaf_force = 0.08,
                              edge_force = 0.25,
                              anchor_force = 0.02,
                              internal_anchor = 0.08,
                              outward_force = 0.01,
                              step = 0.2,
                              cooling = 0.98,
                              tol = 1e-4,
                              ...) {
    initial_layout <- match.arg(initial_layout, c("daylight", "equal_angle", "ape"))
    max_iter <- as.integer(max_iter)
    if (is.na(max_iter) || max_iter < 1L) {
        stop("`max_iter` must be a positive integer.")
    }

    tree_df <- switch(
        initial_layout,
        daylight = suppressMessages(layoutDaylight(model, branch.length, MAX_COUNT = MAX_COUNT)),
        equal_angle = layoutEqualAngle(model, branch.length),
        ape = layoutApe(model, branch.length)
    )

    cache <- .treeAndLeafBuildCache(tree_df)
    if (length(cache$tip_rows) <= 1L) {
        tree_df <- as_tibble(tree_df)
        class(tree_df) <- c("tbl_tree", class(tree_df))
        return(tree_df)
    }

    positions <- as.matrix(tree_df[, c("x", "y")])

    ## The relaxation loop has a compiled version (src/tree_and_leaf.cpp).
    ## Forces here are sums of element-wise contributions, so unlike the
    ## daylight arc accumulation the arithmetic is order-independent; the
    ## R loop remains as reference under the same option switch.
    if (.daylight_use_cpp()) {
        relaxed <- .treeAndLeafRelax_cpp(
            positions, cache,
            list(leaf_force = leaf_force, edge_force = edge_force,
                 anchor_force = anchor_force, internal_anchor = internal_anchor,
                 outward_force = outward_force, step = step,
                 cooling = cooling, tol = tol, max_iter = max_iter)
        )
        tree_df$x <- relaxed[, 1]
        tree_df$y <- relaxed[, 2]
        tree_df <- as_tibble(tree_df)
        class(tree_df) <- c("tbl_tree", class(tree_df))
        return(tree_df)
    }

    initial_positions <- positions
    current_step <- step
    n_nodes <- nrow(tree_df)

    for (iter in seq_len(max_iter)) {
        forces <- matrix(0, nrow = n_nodes, ncol = 2)

        tip_rows <- cache$tip_rows
        tip_weights <- 0.5 + cache$depth_norm[tip_rows]
        n_tips <- length(tip_rows)
        for (i in seq_len(n_tips - 1L)) {
            row_i <- tip_rows[i]
            other_rows <- tip_rows[(i + 1L):n_tips]
            delta <- sweep(positions[other_rows, , drop = FALSE], 2, positions[row_i, ], FUN = "-")
            dist2 <- rowSums(delta^2)
            zero_idx <- dist2 < 1e-12
            if (any(zero_idx)) {
                delta[zero_idx, ] <- .treeAndLeafFallbackDirections(
                    cache$node_ids[row_i],
                    cache$node_ids[other_rows[zero_idx]]
                )
                dist2[zero_idx] <- rowSums(delta[zero_idx, , drop = FALSE]^2)
            }
            dist <- sqrt(dist2)
            magnitude <- leaf_force * (tip_weights[i] + tip_weights[(i + 1L):n_tips]) / (dist2 + 1e-6)
            contribution <- delta / dist * magnitude
            forces[row_i, ] <- forces[row_i, ] - colSums(contribution)
            forces[other_rows, ] <- forces[other_rows, ] + contribution
        }

        for (edge_idx in seq_along(cache$edge_child_rows)) {
            parent_row <- cache$edge_parent_rows[edge_idx]
            child_row <- cache$edge_child_rows[edge_idx]
            delta <- positions[child_row, ] - positions[parent_row, ]
            dist2 <- sum(delta^2)
            if (dist2 < 1e-12) {
                delta <- .treeAndLeafFallbackDirections(
                    cache$node_ids[parent_row],
                    cache$node_ids[child_row]
                )[1, ]
                dist2 <- sum(delta^2)
            }
            dist <- sqrt(dist2)
            magnitude <- edge_force * (dist - cache$target_edge_length[edge_idx])
            contribution <- delta / dist * magnitude
            forces[parent_row, ] <- forces[parent_row, ] + contribution
            forces[child_row, ] <- forces[child_row, ] - contribution
        }

        anchor_weights <- rep(anchor_force, n_nodes)
        anchor_weights[cache$internal_rows] <- internal_anchor
        anchor_weights[cache$root_row] <- 0
        forces <- forces + (initial_positions - positions) * anchor_weights

        if (outward_force != 0) {
            outward_scale <- outward_force * (0.5 + cache$depth_norm)
            forces[cache$tip_rows, ] <- forces[cache$tip_rows, , drop = FALSE] +
                cache$outward_dir[cache$tip_rows, , drop = FALSE] * outward_scale[cache$tip_rows]
        }

        displacement <- forces * (current_step * cache$mobility)
        displacement_norm <- sqrt(rowSums(displacement^2))
        limit_idx <- displacement_norm > current_step & displacement_norm > 0
        if (any(limit_idx)) {
            displacement[limit_idx, ] <- displacement[limit_idx, , drop = FALSE] *
                (current_step / displacement_norm[limit_idx])
        }

        positions <- positions + displacement
        positions[cache$root_row, ] <- initial_positions[cache$root_row, ]

        if (max(displacement_norm, na.rm = TRUE) < tol) {
            break
        }
        current_step <- current_step * cooling
    }

    tree_df$x <- positions[, 1]
    tree_df$y <- positions[, 2]
    tree_df <- as_tibble(tree_df)
    class(tree_df) <- c("tbl_tree", class(tree_df))
    tree_df
}

.treeAndLeafBuildCache <- function(df) {
    node_ids <- as.integer(df$node)
    parent_ids <- as.integer(df$parent)
    n_nodes <- nrow(df)
    max_node <- max(node_ids, na.rm = TRUE)
    row_index <- rep.int(NA_integer_, max_node)
    row_index[node_ids] <- seq_len(n_nodes)

    root_row <- which(is.na(parent_ids) | parent_ids == node_ids)[1]
    if (length(root_row) == 0L || is.na(root_row)) {
        stop("Unable to determine the tree root for `tree_and_leaf` layout.")
    }

    parent_rows <- row_index[parent_ids]
    edge_child_rows <- setdiff(seq_len(n_nodes), root_row)
    edge_parent_rows <- parent_rows[edge_child_rows]

    target_edge_length <- sqrt(
        (df$x[edge_child_rows] - df$x[edge_parent_rows])^2 +
        (df$y[edge_child_rows] - df$y[edge_parent_rows])^2
    )
    target_edge_length[target_edge_length < 1e-8] <- 1e-8

    edge_length_by_row <- numeric(n_nodes)
    edge_length_by_row[edge_child_rows] <- target_edge_length
    depth <- rep(NA_real_, n_nodes)
    depth[root_row] <- 0
    unresolved <- setdiff(seq_len(n_nodes), root_row)
    while (length(unresolved) > 0L) {
        progressed <- FALSE
        for (row in unresolved) {
            parent_row <- parent_rows[row]
            if (!is.na(depth[parent_row])) {
                depth[row] <- depth[parent_row] + edge_length_by_row[row]
                progressed <- TRUE
            }
        }
        unresolved <- which(is.na(depth))
        if (!progressed) {
            depth[is.na(depth)] <- 0
            break
        }
    }
    max_depth <- max(depth, na.rm = TRUE)
    depth_norm <- if (max_depth > 0) depth / max_depth else rep(0, n_nodes)

    outward_dir <- cbind(df$x - df$x[root_row], df$y - df$y[root_row])
    dir_norm <- sqrt(rowSums(outward_dir^2))
    fallback_rows <- dir_norm < 1e-12
    if (any(fallback_rows)) {
        fallback <- .treeAndLeafFallbackDirections(node_ids[fallback_rows], node_ids[fallback_rows] + 1L)
        outward_dir[fallback_rows, ] <- fallback
        dir_norm[fallback_rows] <- sqrt(rowSums(fallback^2))
    }
    outward_dir <- outward_dir / dir_norm

    tip_rows <- which(df$isTip)
    internal_rows <- which(!df$isTip)
    mobility <- rep(0.18, n_nodes)
    mobility[tip_rows] <- 0.75 + 0.25 * depth_norm[tip_rows]
    mobility[root_row] <- 0

    list(
        node_ids = node_ids,
        root_row = root_row,
        tip_rows = tip_rows,
        internal_rows = internal_rows,
        edge_parent_rows = edge_parent_rows,
        edge_child_rows = edge_child_rows,
        target_edge_length = target_edge_length,
        depth_norm = depth_norm,
        outward_dir = outward_dir,
        mobility = mobility
    )
}

.treeAndLeafFallbackDirections <- function(node_ids, other_ids) {
    angles <- (((as.numeric(node_ids) + as.numeric(other_ids)) * 0.618033988749895) %% 1) * 2 * pi
    cbind(cos(angles), sin(angles))
}

##' Apply the daylight alorithm to adjust the spacing between the subtrees and tips of the
##' specified node.
##'
##' @title applyLayoutDaylight
##' @param df tree data.frame
##' @param node_id is id of the node from which daylight is measured to the other subtrees.
##' @param cache optional precomputed daylight cache.
##' @return list with tree data.frame with updated layout using daylight algorithm and max_change angle.
##' @importFrom rlang .data
##
##
## ```
## for node in nodes {
##   if node is a leaf {
##     next
##   }
##
##   subtrees = get subtrees of node
##
##   for i-th subtree in subtrees {
##     [end, start] = get left and right angles of tree from node id.
##     angle_list[i, 'left'] = end
##     angle_list[i, 'beta'] = start - end  # subtree arc angle
##     angle_list[i, 'index'] = i-th # index of subtree/leaf
##   }
##
##   sort angle_list by 'left' column in ascending order.
##
##   D = 360 - sum( angle_list['beta'] ) # total daylight angle
##   d = D / |subtrees| # equal daylight angle.
##
##   new_L = left angle of first subtree.
##
##   for n-th row in angle_list{
##     # Calculate angle to rotate subtree/leaf to create correct daylight angle.
##     new_left_angle = new_left_angle + d + angle_list[n, 'beta']
##     Calculate the difference between the old and new left angles.
##     adjust_angle = new_left_angle - angle_list[n, 'left']
##
##     index = angle_list['index']
##     rotate subtree[index] wrt n-th node by adjust_angle
##     }
##   }
## }
## ```
##' @noRd
applyLayoutDaylight <- function(df, node_id, cache = NULL){
  if (is.null(cache)) {
    cache <- .daylightBuildCache(df)
  }

  subtrees <- cache$subtrees[[node_id]]
  if (is.null(subtrees) || length(subtrees) <= 2L) {
    return(list(tree = df, max_change = 0.0))
  }

  n_subtrees <- length(subtrees)
  angle_left <- numeric(n_subtrees)
  angle_beta <- numeric(n_subtrees)

  ## One reusable theta buffer per layout call. Filled per subtree rather than
  ## per origin, because rotating the earlier subtrees moves the later ones and
  ## changes their angles. The buffer is reused, so the allocation that used to
  ## dominate (2996 calls x 16 KB at 1000 tips) is gone.
  if (.daylight_use_cpp()) {
    if (is.null(cache$df_x_num)) {
      cache$row_index_int <- as.integer(cache$row_index)
      cache$has_children_lgl <- as.logical(cache$has_children)
      cache$theta_buffer <- numeric(nrow(df))
      cache$df_x_num <- TRUE   # marks the fields above as initialised
    }
  }

  for (i in seq_len(n_subtrees)) {
    arc <- .daylightGetArcAngles(df, node_id, subtrees[[i]], cache)
    angle_left[i] <- arc[['left']]
    beta <- arc[['left']] - arc[['right']]
    if (beta < 0) {
      beta <- beta + 2
    }
    angle_beta[i] <- beta
  }

  order_idx <- order(angle_left)
  angle_left <- angle_left[order_idx]
  angle_beta <- angle_beta[order_idx]

  total_daylight <- 2 - sum(angle_beta)
  d <- total_daylight / n_subtrees
  new_left_angle <- angle_left[1]

  max_change <- 0.0
  for (i in 2:n_subtrees) {
    new_left_angle <- new_left_angle + d + angle_beta[i]
    adjust_angle <- new_left_angle - angle_left[i]
    max_change <- max(max_change, abs(adjust_angle))
    subtree_nodes <- subtrees[[order_idx[i]]]$subtree
    df <- .daylightRotatePoints(df, node_id, subtree_nodes, adjust_angle, cache)
  }

  list(tree = df, max_change = max_change)
}

.daylightBuildCache <- function(df) {
  node_ids <- as.integer(df$node)
  parent_ids <- as.integer(df$parent)
  max_node <- max(node_ids, na.rm = TRUE)

  row_index <- rep.int(NA_integer_, max_node)
  row_index[node_ids] <- seq_along(node_ids)

  children <- vector('list', max_node)
  has_parent <- !is.na(parent_ids) & parent_ids != node_ids
  if (any(has_parent)) {
    split_children <- split(node_ids[has_parent], parent_ids[has_parent])
    for (i in seq_along(split_children)) {
      parent_id <- as.integer(names(split_children)[i])
      children[[parent_id]] <- unname(as.integer(split_children[[i]]))
    }
  }
  empty_children <- lengths(children) == 0L
  children[empty_children] <- list(integer())

  has_children <- lengths(children) > 0L
  is_tip <- rep.int(FALSE, max_node)
  is_tip[node_ids] <- !node_ids %in% parent_ids[has_parent]

  root_rows <- which(is.na(df$parent) | df$parent == df$node)
  root_node <- if (length(root_rows) > 0L) node_ids[root_rows[1]] else NA_integer_

  bfs_nodes <- integer()
  if (!is.na(root_node)) {
    bfs_nodes <- root_node
    i <- 1L
    while (i <= length(bfs_nodes)) {
      child_ids <- children[[bfs_nodes[i]] ]
      i <- i + 1L
      if (length(child_ids) > 0L) {
        bfs_nodes <- c(bfs_nodes, child_ids)
      }
    }
  }
  internal_nodes <- bfs_nodes[has_children[bfs_nodes]]

  subtrees <- vector('list', max_node)
  for (node_id in internal_nodes) {
    child_ids <- children[[node_id]]
    node_subtrees <- lapply(child_ids, function(child_id) {
      list(node = child_id, subtree = getSubtree.df(df, child_id))
    })

    covered_nodes <- unlist(lapply(node_subtrees, `[[`, 'subtree'), use.names = FALSE)
    remaining_nodes <- setdiff(node_ids, covered_nodes)
    parent_id <- parent_ids[row_index[node_id]]
    if (!is.na(parent_id) && parent_id != node_id && length(remaining_nodes) > 0L) {
      node_subtrees[[length(node_subtrees) + 1L]] <- list(node = parent_id, subtree = remaining_nodes)
    }
    subtrees[[node_id]] <- node_subtrees
  }

  list(
    row_index = row_index,
    children = children,
    has_children = has_children,
    is_tip = is_tip,
    subtrees = subtrees,
    root_node = root_node,
    internal_nodes = internal_nodes
  )
}

.daylightGetArcAngles <- function(df, origin_id, subtree, cache) {
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

  ## The accumulation loop is the hot spot; it also has a compiled version
  ## that is bit-identical (see src/daylight_arc.cpp and the parity test).
  ## theta_left/theta_right are already resolved at this point, so the switch
  ## covers only the loop.
  if (.daylight_use_cpp()) {
    ## theta_by_id is prepared once per origin by applyLayoutDaylight() (see
    ## .daylight_theta_cache), so the subtrees of the same origin share it
    ## instead of each allocating a full-length scratch vector.
    if (!is.null(cache$theta_buffer)) {
      theta <- cache$theta_buffer
      daylight_fill_theta_cpp(
        df_x = as.numeric(df$x), df_y = as.numeric(df$y),
        row_index = cache$row_index_int,
        origin_id = as.integer(origin_id),
        node_ids = as.integer(subtree$subtree),
        theta_by_id = theta
      )
      return(daylight_arc_angles_cached_cpp(
        df_x = as.numeric(df$x), df_y = as.numeric(df$y),
        row_index = cache$row_index_int,
        origin_id = as.integer(origin_id),
        subtree_node_ids = as.integer(subtree$subtree),
        children = cache$children,
        has_children = cache$has_children_lgl,
        theta_left_start = as.numeric(theta_left),
        theta_right_start = as.numeric(theta_right),
        theta_by_id = theta
      ))
    }
    return(.daylight_arc_cpp(df, origin_id, subtree, cache, theta_left, theta_right))
  }

  arc <- c(left = theta_left, right = theta_right)
  parent_nodes <- subtree_node_ids[subtree_node_ids != origin_id]
  parent_nodes <- parent_nodes[cache$has_children[parent_nodes]]

  ## Angle of every node the loop below looks at, in one vectorised atan2().
  ##
  ## atan2() is a pure function, so this is bit-identical to calling
  ## getNodeAngle.vector() per node, which is what dominated the profile: the
  ## inner loop ran once per descendant per subtree, per round.
  ##
  ## Only the *angles* are hoisted. The loop itself keeps its original order:
  ## `arc` is updated in place and every comparison reads the value written by
  ## the previous iteration, so reordering it or replacing the accumulation
  ## with max()/min() would change the result.
  children_of <- cache$children
  looked_at <- unique(c(parent_nodes,
                        unlist(children_of[parent_nodes], use.names = FALSE)))
  theta_by_id <- rep.int(NA_real_, length(row_index))
  if (length(looked_at) > 0L) {
    rows <- row_index[looked_at]
    theta_by_id[looked_at] <- getNodeAngle.vector(x_origin, y_origin,
                                                  df_x[rows], df_y[rows])
  }

  for (parent_id in parent_nodes) {
    theta_parent <- theta_by_id[parent_id]
    child_ids <- children_of[[parent_id]]
    child_ids <- child_ids[child_ids != origin_id]
    if (length(child_ids) == 0L) {
      next
    }

    for (child_id in child_ids) {
      theta_child <- theta_by_id[child_id]
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

##'
##' One call instead of three. Returns the left and right boundary angles of
##' the clade subtree, the largest distance from `node` to any node inside
##' the clade, and the largest distance from the root of the whole tree to
##' any node -- the two radii the clade label / hilight code needs.
##'
##' @param data a node table with columns `node`, `parent`, `x`, `y`.
##' @param node the root node id of the clade.
##' @return a list with `left`, `right`, `radius` and `root_radius`.
##'
##' @examples
##' tr <- ape::rtree(10)
##' d <- tree_layout(tr, layout = "circular")
##' clade_geometry(d, node = 12)
##' @export
clade_geometry <- function(data, node) {
    subtree <- list(subtree = getSubtree.df(data, node), node = node)
    arc <- getTreeArcAngles(data, node, subtree)
    in_clade <- data$node %in% subtree$subtree
    list(
        left        = as.numeric(arc["left"]),
        right       = as.numeric(arc["right"]),
        radius      = max(getNodeEuclDistances(data[in_clade, ], node)),
        root_radius = max(getNodeEuclDistances(data, .getRoot.df(data)))
    )
}

##' Find the right (clockwise rotation, angle from +ve x-axis to furthest subtree nodes) and
##' left (anti-clockwise angle from +ve x-axis to subtree) Returning arc angle in `[0, 2]` (0 to 360) domain.
##'
##' @title getTreeArcAngles
##' @param df tree data.frame
##' @param origin_id node id from which to calculate left and right hand angles of subtree.
##' @param subtree named list of root id of subtree (node) and list of node ids for given subtree (subtree).
##' @return named list with right and left angles in range `[0,2]` i.e 1 = 180 degrees, 1.5 = 270 degrees.
##' @noRd

getTreeArcAngles <- function(df, origin_id, subtree) {
  .daylightGetArcAngles(df, origin_id, subtree, .daylightBuildCache(df))
}

.daylightRotatePoints <- function(df, pivot_node, nodes, angle, cache) {
  node_rows <- cache$row_index[nodes]
  node_rows <- node_rows[!is.na(node_rows)]
  if (length(node_rows) == 0L) {
    return(df)
  }

  pivot_row <- cache$row_index[pivot_node]
  pivot_x <- df$x[pivot_row]
  pivot_y <- df$y[pivot_row]
  cospitheta <- cospi(angle)
  sinpitheta <- sinpi(angle)
  delta_x <- df$x[node_rows] - pivot_x
  delta_y <- df$y[node_rows] - pivot_y

  df$x[node_rows] <- cospitheta * delta_x - sinpitheta * delta_y + pivot_x
  df$y[node_rows] <- sinpitheta * delta_x + cospitheta * delta_y + pivot_y

  tip_nodes <- nodes[cache$is_tip[nodes]]
  if (length(tip_nodes) > 0L) {
    tip_rows <- cache$row_index[tip_nodes]
    parent_rows <- cache$row_index[df$parent[tip_rows]]
    theta <- getNodeAngle.vector(df$x[parent_rows], df$y[parent_rows], df$x[tip_rows], df$y[tip_rows])
    df$angle[tip_rows] <- 180 * ifelse(theta < 0, 2 + theta, theta)
  }

  df
}
##'
##' @title getNodeAngle.df
##' @param df tree data.frame
##' @param origin_node_id origin node id number
##' @param node_id end node id number
##' @return angle in range `[-1, 1]`, i.e. degrees/180, radians/pi
##' @noRd
getNodeAngle.df <- function(df, origin_node_id, node_id){
  if (origin_node_id != node_id) {
    df_x = df$x
    df_y = df$y
    atan2(df_y[node_id] - df_y[origin_node_id], df_x[node_id] - df_x[origin_node_id]) / pi
  }else{
    NA
  }
}

getNodeAngle.vector <- function(x_origin, y_origin, x, y) {
  atan2(y - y_origin, x - x_origin) / pi
}

euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

## Get the distances from the node to all other nodes in data.frame (including itself if in df)
##
## Vectorised.  The previous per-row loop re-subset the whole `df` inside
## every iteration and called `euc.dist()` once per node, which is quadratic
## (~5.5 s on a 4000-tip tree).  The arithmetic below is identical.
##' Euclidean distances from a node
##'
##' @keywords internal
##' @noRd
getNodeEuclDistances <- function(df, node){
  # https://stackoverflow.com/questions/24746892/how-to-calculate-euclidian-distance-between-two-points-defined-by-matrix-contain#24747155
  i <- match(node, df$node)
  if (is.na(i)) {
      stop("`node` is not found in `df$node`.", call. = FALSE)
  }
  sqrt((df$x - df$x[i])^2 + (df$y - df$y[i])^2)
}


##' Get all children of node from tree, including start_node.
##'
##' @title getSubtree
##' @param tree ape phylo tree object
##' @param node is the tree node id from which the tree is derived.
##' @return list of all child node id's from starting node.
##' @noRd
getSubtree <- function(tree, node){

  ## subtree <- c(node)
  ## i <- 1
  ## while( i <= length(subtree)){
  ##   subtree <- c(subtree, tidytree::child(tree, subtree[i]))
  ##   # remove any '0' root nodes
  ##   subtree <- subtree[subtree != 0]
  ##   i <- i + 1
  ## }
    ## return(subtree)
    offspring(tree, node, self_include = TRUE)
}

##' Get all children of node from df tree using breath-first.
##'
##' @title getSubtree.df
##' @param df tree data.frame
##' @param node id of starting node.
##' @return list of all child node id's from starting node.
##' @noRd
getSubtree.df <- function(df, node){
  ## subtree <- node[node != 0]
  ## i <- 1
  ## while( i <= length(subtree)){
  ##     ## subtree <- c(subtree, getChild.df(df, subtree[i]))
  ##     subtree <- c(subtree, tidytree::child(df, subtree[i])$node)
  ##   i <- i + 1
  ## }
    ## subtree
    #tidytree:::offspring.tbl_tree(df, node, self_include = TRUE)$node
    offspring(df, node, self_include = TRUE)$node
}

##' Get all subtrees of specified node. This includes all ancestors and relatives of node and
##' return named list of subtrees.
##'
##' @title getSubtreeUnrooted
##' @param tree ape phylo tree object
##' @param node is the tree node id from which the subtrees are derived.
##' @return named list of subtrees with the root id of subtree and list of node id's making up subtree.
##' @noRd
getSubtreeUnrooted <- function(tree, node){
  # if node leaf, return nothing.
  if( isTip(tree, node) ){
    # return NA
    return(NA)
  }

  subtrees <- list()

  # get subtree for each child node.
    ## children_ids <- getChild(tree, node)
    children_ids <- child(tree, node)

  remaining_nodes <- getNodes_by_postorder(tree)
  # Remove current node from remaining_nodes list.
  remaining_nodes <- setdiff(remaining_nodes, node)


  for( child in children_ids ){
    # Append subtree nodes to list if not 0 (root).
    subtree <- getSubtree(tree, child)
    subtrees[[length(subtrees)+1]] <- list( node = child, subtree = subtree)
    # remove subtree nodes from remaining nodes.
    remaining_nodes <- setdiff(remaining_nodes, as.integer(unlist(subtrees[[length(subtrees)]]['subtree']) ))
  }

  # The remaining nodes that are not found in the child subtrees are the remaining subtree nodes.
  # ie, parent node and all other nodes. We don't care how they are connect, just their ids.
  parent_id <- parent(tree, node)
  # If node is not root, add remainder of tree nodes as subtree.
  if( parent_id != 0 & length(remaining_nodes) >= 1){
    subtrees[[length(subtrees)+1]] <- list( node = parent_id, subtree = remaining_nodes)
  }

  return(subtrees)
}



##' @keywords internal
##' @noRd
.getRoot.df <- function(df, node){

  root <- which(is.na(df$parent))
  # Check if root was found.
  if(length(root) == 0){
      ## Alternatively, root can self reference, eg node = 10, parent = 10
      root <- df$node[df$parent == df$node]
      ## root <- unlist(apply(df, 1, function(x){ if(x['node'] == x['parent']){ x['node'] } }))
  }
  return(root)
}




##' Node ids in post-order
##'
##' @param tree a `phylo` object.
##' @return an integer vector of node ids.
##'
##' @importFrom ape reorder.phylo
##' @keywords internal
##' @noRd
getNodes_by_postorder <- function(tree) {
  tree <- reorder.phylo(tree, "postorder")
    unique(rev(as.vector(t(tree$edge[,c(2,1)]))))
}

##' x coordinates from an edge list
##'
##' Computes the x coordinate of every node by walking the edge list from a
##' starting root position. Used by code that recomputes horizontal placement
##' after changing the tree order.
##'
##' @param x numeric vector indexed by node id; the first entry is overwritten
##'   with `start`.
##' @param root node id to start from.
##' @param parent parent node ids of the edge list.
##' @param child child node ids of the edge list.
##' @param len branch length of every edge.
##' @param start x coordinate assigned to the root, default 0.
##' @param rev logical; if TRUE the tree is walked in reverse, so coordinates
##'   decrease from the root.
##' @return the `x` vector with one coordinate per node.
##' @examples
##' tr <- ape::rtree(6)
##' x <- numeric(max(tr$edge))
##' getXcoord2(x, root = ape::Ntip(tr) + 1L,
##'            parent = tr$edge[, 1], child = tr$edge[, 2],
##'            len = tr$edge.length)
##' @export
getXcoord2 <- function(x, root, parent, child, len, start=0, rev=FALSE) {
    direction <- if (isTRUE(rev)) -1 else 1
    ignore_negative_edge <- getOption("ignore.negative.edge", default=FALSE)

    if (any(len < 0) && !ignore_negative_edge) {
        warning("The tree contained negative ",
                ifelse(sum(len < 0) > 1, "edge lengths", "edge length"), ". ",
                "If you want to ignore the ",
                ifelse(sum(len < 0) > 1, "edges", "edge"), ", you can set ",
                "options(ignore.negative.edge=TRUE), then re-run ggtree.")
    }

    step <- if (ignore_negative_edge) {
        len * direction * sign(len)
    } else {
        len * direction
    }

    ## One top-down pass.
    ##
    ## The previous loop advanced one level per iteration (`parent %in% current`),
    ## so a ladder-like tree -- depth ~ n -- cost O(n^2) and `slanted` took 4.8 s
    ## at 2000 tips. This only has vectors, not a tree object, so the order is
    ## derived from the edge list itself: repeatedly take nodes whose parent
    ## already has a value, which is one pass when walked breadth-first.
    x[root] <- start
    x[-root] <- NA

    ## children grouped by parent
    max_id <- max(child, parent, root)
    cl <- split(child, parent)
    kids_of <- vector("list", max_id)
    kids_of[as.integer(names(cl))] <- cl
    cl_len <- split(step, parent)
    len_of <- vector("list", max_id)
    len_of[as.integer(names(cl_len))] <- cl_len

    frontier <- root
    while (length(frontier) > 0L) {
        nxt <- vector("list", length(frontier))
        for (i in seq_along(frontier)) {
            v <- frontier[i]
            kids <- kids_of[[v]]
            if (length(kids) > 0L) {
                x[kids] <- x[v] + len_of[[v]]
                nxt[[i]] <- kids
            }
        }
        frontier <- unlist(nxt, use.names = FALSE)
    }

    return(x)
}






getXcoord_no_length <- function(tr) {
    edge <- tr$edge
    parent <- edge[,1]
    child <- edge[,2]

    N <- getNodeNum(tr)
    ntip <- Ntip(tr)
    x <- numeric(N)
    ## every internal node sits one unit left of its leftmost child, so the
    ## tips (which have none) stay at 0 until the shift at the end

    cl <- split(child, parent)
    child_list <- list()
    child_list[as.numeric(names(cl))] <- cl

    ## One upward pass, same reason as getYcoord(): the previous loop advanced
    ## one level per iteration and rebuilt table()/factor() structures each
    ## time, which on a ladder-like tree (depth ~ n) made the cost grow faster
    ## than the node count.
    for (v in order(ape::node.depth(tr))) {
        ch <- child_list[[v]]
        if (length(ch) > 0L) {
            x[v] <- min(x[ch], na.rm = TRUE) - 1
        }
    }

    x <- x - min(x, na.rm = TRUE)
    return(x)
}




getXcoord <- function(tr) {
    edge <- tr$edge
    parent <- edge[,1]
    child <- edge[,2]
    root <- getRoot(tr)

    len <- tr$edge.length

    N <- getNodeNum(tr)
    x <- numeric(N)
    x <- getXcoord2(x, root, parent, child, len)
    return(x)
}




## scale the branch (the line plotted) to the actual value of edge length
## but it seems not the good idea as if we want to add x-axis (e.g. time-scaled tree)
## then the x-value is not corresponding to edge length as in rectangular layout
## getXYcoord_slanted <- function(tr) {
##     edge <- tr$edge
##     parent <- edge[,1]
##     child <- edge[,2]
##     root <- getRoot(tr)

##     N <- getNodeNum(tr)
##     len <- tr$edge.length
##     y <- getYcoord(tr, step=min(len)/2)
##     len <- sqrt(len^2 - (y[parent]-y[child])^2)
##     x <- numeric(N)
##     x <- getXcoord2(x, root, parent, child, len)
##     res <- data.frame(x=x, y=y)
##     return(res)
## }

## Tidy tree layout, from Ploeg (2014) Softw. Pract. Exper. 44:1467-1484.
##
## NOTE on multifurcations: Ploeg's `distributeExtra` / `addChildSpacing`
## steps, which spread the extra separation over the *middle* children of a
## node with 3 or more children, are not implemented here.  Binary trees are
## unaffected; polytomous trees are still laid out without overlaps, but the
## middle children of a polytomy are not optimally spaced.
##
## NOTE on tip separation: this layout allows two tips to share the same `y`
## when they sit at different depths.  The drawing is still valid (no branch
## passes through a node) but `y` is no longer a unique key per tip, so
## y-based annotations such as `gheatmap()` will overlap.

tidy_positionRoot <- function(node, ch, nch, coords.env) {
  coords.env$prelim[node] <- (coords.env$prelim[ch[1]] + coords.env$mod[ch[1]] +
                                coords.env$prelim[ch[nch]] + coords.env$mod[ch[nch]] +
                                coords.env$w[ch[nch]]) / 2 -
    coords.env$w[node] / 2
}

tidy_moveSubtree <- function(ch, i, dist, coords.env) {
  coords.env$mod[ch[i]]  <- coords.env$mod[ch[i]]  + dist
  coords.env$msel[ch[i]] <- coords.env$msel[ch[i]] + dist
  coords.env$mser[ch[i]] <- coords.env$mser[ch[i]] + dist
}

tidy_nextLeftContour <- function(node, coords.env) {
  ch <- coords.env$children[[node]]
  if (length(ch) == 0) coords.env$tl[node] else ch[1L]
}

tidy_nextRightContour <- function(node, coords.env) {
  ch <- coords.env$children[[node]]
  if (length(ch) == 0) coords.env$tr[node] else ch[length(ch)]
}

tidy_setLeftThread <- function(ch, i, cl, modsumcl, coords.env) {
  li <- coords.env$el[ch[1]]
  coords.env$tl[li] <- cl
  diff <- modsumcl - coords.env$mod[cl] - coords.env$msel[ch[1]]
  coords.env$mod[li] <- coords.env$mod[li] + diff
  coords.env$prelim[li] <- coords.env$prelim[li] - diff
  coords.env$el[ch[1]] <- coords.env$el[ch[i]]
  coords.env$msel[ch[1]] <- coords.env$msel[ch[i]]
}

tidy_setRightThread <- function(ch, i, sr, modsumsr, coords.env) {
  ri <- coords.env$er[ch[i]]
  coords.env$tr[ri] <- sr
  diff <- modsumsr - coords.env$mod[sr] - coords.env$mser[ch[i]]
  coords.env$mod[ri] <- coords.env$mod[ri] + diff
  coords.env$prelim[ri] <- coords.env$prelim[ri] - diff
  coords.env$er[ch[i]] <- coords.env$er[ch[i-1]]
  coords.env$mser[ch[i]] <- coords.env$mser[ch[i-1]]
}

tidy_separate <- function(ch, i, coords.env) {
  sr <- ch[i - 1];
  mssr <- coords.env$mod[sr]
  cl <- ch[i];
  mscl <- coords.env$mod[cl]

  while (!is.na(sr) && !is.na(cl)) {
    dist <- mssr + coords.env$prelim[sr] + coords.env$w[sr] - (mscl + coords.env$prelim[cl])

    if (dist > 0) {
      mscl <- mscl + dist

      tidy_moveSubtree(ch, i, dist, coords.env)
    }

    sy <- coords.env$y[sr]
    cy <- coords.env$y[cl]

    if (sy <= cy) {
      sr <- tidy_nextRightContour(sr, coords.env)

      if (!is.na(sr)) mssr <- mssr + coords.env$mod[sr]
    }

    if (sy >= cy) {
      cl <- tidy_nextLeftContour(cl, coords.env)

      if (!is.na(cl)) mscl <- mscl + coords.env$mod[cl]
    }
  }

  if (is.na(sr) && !is.na(cl)) {
    tidy_setLeftThread(ch, i, cl, mscl, coords.env)
  } else if (!is.na(sr) && is.na(cl)) {
    tidy_setRightThread(ch, i, sr, mssr, coords.env)
  }
}

## Both walks are iterative (explicit stack).  The recursive versions
## overflowed R's node stack on trees deeper than ~1500 nodes; this also
## removes the per-node rescan of the edge matrix.
tidy_first_walk <- function(root, coords.env) {
  children <- coords.env$children
  N <- length(children)

  ## pre-order traversal; reversing it yields children before parents
  ord <- integer(N)
  stk <- integer(N)
  k <- 0L
  sp <- 1L
  stk[1L] <- root
  while (sp > 0L) {
    v <- stk[sp]
    sp <- sp - 1L
    k <- k + 1L
    ord[k] <- v
    ch <- children[[v]]
    nc <- length(ch)
    if (nc > 0L) {
      stk[(sp + 1L):(sp + nc)] <- rev(ch)
      sp <- sp + nc
    }
  }

  for (j in k:1L) {
    node <- ord[j]
    ch <- children[[node]]
    nch <- length(ch)

    if (nch == 0L) {                      ## setExtremes
      coords.env$el[node]   <- node
      coords.env$er[node]   <- node
      coords.env$msel[node] <- 0
      coords.env$mser[node] <- 0
      next
    }

    for (i in seq_along(ch)[-1]) {
      tidy_separate(ch, i, coords.env)
    }

    tidy_positionRoot(node, ch, nch, coords.env)

    coords.env$el[node]   <- coords.env$el[ch[1]]
    coords.env$msel[node] <- coords.env$msel[ch[1]]
    coords.env$er[node]   <- coords.env$er[ch[nch]]
    coords.env$mser[node] <- coords.env$mser[ch[nch]]
  }
}

tidy_second_walk <- function(root, coords.env) {
  children <- coords.env$children
  N <- length(children)

  ns <- integer(N)
  ms <- numeric(N)
  sp <- 1L
  ns[1L] <- root
  ms[1L] <- 0
  while (sp > 0L) {
    node <- ns[sp]
    modsum <- ms[sp]
    sp <- sp - 1L

    modsum <- modsum + coords.env$mod[node]
    coords.env$x[node] <- coords.env$prelim[node] + modsum

    ch <- children[[node]]
    nc <- length(ch)
    if (nc > 0L) {
      ns[(sp + 1L):(sp + nc)] <- rev(ch)
      ms[(sp + 1L):(sp + nc)] <- modsum
      sp <- sp + nc
    }
  }
}

##' Tidy (non-layered) tree layout
##'
##' Computes the vertical (`y`) coordinate by the linear-time algorithm of
##' Ploeg (2014).  The horizontal (`x`) coordinate is taken as given, so it
##' plays the role of the "level" in the original algorithm -- hence `x` and
##' `y` are swapped in the environment below.
##' @param tr phylo object
##' @param x horizontal coordinate (depth) of every node
##' @param step minimum separation between adjacent nodes
##' @param extra.tip.padding amount, relative to `max(x)`, by which tips are
##' pushed away from the root in the level comparison.  Larger values keep
##' tips further apart, at the cost of a taller drawing.
##' @return vertical coordinate of every node
##' @references Ploeg, A. (2014) Drawing non-layered tidy trees in linear time.
##' \emph{Software: Practice and Experience} 44(12):1467-1484.
##' doi:10.1002/spe.2213
##' @noRd
##' @importFrom rlang env
getYcoord_tidy <- function(tr, x, step = 1, extra.tip.padding = 0.01) {
  Ntip <- length(tr$tip.label)
  N <- getNodeNum(tr)
  root <- getRoot(tr)

  ## children indexed by node number, in edge order, so the walks never have
  ## to rescan the edge matrix
  children <- split(tr$edge[, 2], factor(tr$edge[, 1], levels = seq_len(N)))

  # environment/output for the walks
  # x and y are swapped
  coords.env <- env(
    children = children,
    w = rep(step, N),
    y = x + c(rep(max(x) * extra.tip.padding, Ntip), rep(0, N - Ntip)),
    # outputs
    prelim = rep(0, N), mod = rep(0, N),
    x = rep(NA, N),
    tl = rep(NA, N), tr = rep(NA, N),
    el = rep(NA, N), er = rep(NA, N),
    msel = rep(NA, N), mser = rep(NA, N)
  )

  tidy_first_walk(root, coords.env)
  tidy_second_walk(root, coords.env)

  # x and y are swapped
  return(coords.env$x)
}

## @importFrom magrittr %>%
##' y coordinates for tips
##'
##' Places tips on evenly spaced y positions (or in the order given by
##' `tip.order`) and fills the internal nodes bottom-up with the mean of their
##' children.
##'
##' @param tr a `phylo` object.
##' @param step spacing between adjacent tips, default 1.
##' @param tip.order tip labels in the order they should appear bottom to top;
##'   `NULL` keeps the existing tip order.
##' @return numeric vector of y coordinates, indexed like the nodes of `tr`.
##' @examples
##' tr <- ape::rtree(6)
##' getYcoord(tr)
##' getYcoord(tr, tip.order = rev(tr$tip.label))
##' @export
getYcoord <- function(tr, step=1, tip.order = NULL) {
    Ntip <- length(tr[["tip.label"]])
    N <- getNodeNum(tr)

    edge <- tr[["edge"]]
    parent <- edge[,1]
    child <- edge[,2]

    cl <- split(child, parent)
    child_list <- list()
    child_list[as.numeric(names(cl))] <- cl

    y <- numeric(N)
    if (is.null(tip.order)) {
        tip.idx <- child[child <= Ntip]
        y[tip.idx] <- 1:Ntip * step
    } else {
        tip.idx <- 1:Ntip
        y[tip.idx] <- match(tr$tip.label, tip.order) * step
    }
    y[-tip.idx] <- NA

    ## One upward pass.
    ##
    ## This used to walk the tree one level per iteration, rescanning the
    ## pending tips with `%in%` each time. A ladder-like tree has depth ~ n, so
    ## the cost grew faster than the number of nodes: 0.16 s at 200 tips,
    ## 2.8 s at 800, 17.7 s at 1600, against 0.07 s at 1600 for a balanced
    ## tree. Every layout that is not unrooted, tidy or external goes through
    ## here.
    ##
    ## ape::node.depth() counts edges from the tips, so increasing order visits
    ## every child before its parent and one pass is enough.
    for (v in order(ape::node.depth(tr))) {
        ch <- child_list[[v]]
        if (length(ch) > 0L) {
            y[v] <- mean(y[ch], na.rm = TRUE)
        }
    }

    return(y)
}



##' y coordinates scaled by a number
##'
##' @keywords internal
##' @noRd
getYcoord_scale_numeric <- function(tr, df, yscale, ...) {
    df <- .assign_parent_status(tr, df, yscale)
    df <- .assign_child_status(tr, df, yscale)

    y <- df[, yscale]

    if (anyNA(y)) {
        warning("NA found in y scale mapping, all were setting to 0")
        y[is.na(y)] <- 0
    }

    return(y)
}


.assign_parent_status <- function(tr, df, variable) {
    yy <- df[[variable]]
    na.idx <- which(is.na(yy))
    if (length(na.idx) > 0) {
        tree <- get.tree(tr)
        nodes <- getNodes_by_postorder(tree)
        ## children once, instead of child(tree, curNode) inside the loop:
        ## that rescanned the whole table on every node, which is O(n^2) on a
        ## ladder-like tree
        kids_of <- .children_index(tree)
        for (curNode in nodes) {
            children <- kids_of[[curNode]]
            if (length(children) == 0) {
                next
            }
            idx <- which(is.na(yy[children]))
            if (length(idx) > 0) {
                yy[children[idx]] <- yy[curNode]
            }
        }
    }
    df[, variable] <- yy
    return(df)
}


.assign_child_status <- function(tr, df, variable, yscale_mapping=NULL) {
    yy <- df[[variable]]
    if (!is.null(yscale_mapping)) {
        yy <- yscale_mapping[yy]
    }

    na.idx <- which(is.na(yy))
    if (length(na.idx) > 0) {
        tree <- get.tree(tr)
        nodes <- rev(getNodes_by_postorder(tree))
        ## both lookups built once; the loop body used to call parent() and
        ## child(), each of which scanned the whole table
        kids_of <- .children_index(tree)
        parent_of <- .parent_index(tree)
        for (curNode in nodes) {
            parent <- parent_of[curNode]
            if (parent == 0) { ## already reach root
                next
            }
            idx <- which(is.na(yy[parent]))
            if (length(idx) > 0) {
                child <- kids_of[[parent]]
                yy[parent[idx]] <- mean(yy[child], na.rm=TRUE)
            }
        }
    }
    df[, variable] <- yy
    return(df)
}


##' y coordinates scaled by a category
##'
##' @keywords internal
##' @noRd
getYcoord_scale_category <- function(tr, df, yscale, yscale_mapping=NULL, ...) {
    if (is.null(yscale_mapping)) {
        stop("yscale is category variable, user should provide yscale_mapping,
             which is a named vector, to convert yscale to numberical values...")
    }
    if (! is(yscale_mapping, "numeric") ||
        is.null(names(yscale_mapping))) {
        stop("yscale_mapping should be a named numeric vector...")
    }

    if (yscale == "label") {
        yy <- df[[yscale]]
        ii <- which(is.na(yy))
        if (length(ii)) {
            ## df[ii, yscale] <- df[ii, "node"]
            df[[yscale]][ii] <- as.character(df[['node']][ii])
        }
    }

    ## assign to parent status is more prefer...
    df <- .assign_parent_status(tr, df, yscale)
    df <- .assign_child_status(tr, df, yscale, yscale_mapping)

    y <- df[[yscale]]

    if (anyNA(y)) {
        warning("NA found in y scale mapping, all were setting to 0")
        y[is.na(y)] <- 0
    }
    return(y)
}


add_angle_slanted <- function(res) {
    x <- res[["x"]]
    y <- res[["y"]]
    dy <- (y - y[match(res$parent, res$node)]) / diff(range(y))
    dx <- (x - x[match(res$parent, res$node)]) / diff(range(x))
    theta <- atan(dy/dx)
    theta[is.na(theta)] <- 0 ## root node
    res$angle <- theta/pi * 180

    branch.y <- (y[match(res$parent, res$node)] + y)/2
    idx <- is.na(branch.y)
    branch.y[idx] <- y[idx]
    res[, "branch.y"] <- branch.y
    return(res)
}


##' Branch mid-point coordinates
##'
##' Adds `branch` (the mid-point x), and for the unrooted layouts also
##' `branch.x` and `branch.y`, which place labels on the middle of a branch
##' even when the branch is not axis-aligned.
##'
##' @param res a node table with `node`, `parent`, `x` and `y` columns.
##' @param layout the layout name the coordinates came from; `branch.x` and
##'   `branch.y` are only computed for the unrooted layouts.
##' @return the node table with the branch columns added.
##' @examples
##' tr <- ape::rtree(6)
##' d <- tree_layout(tr, layout = "equal_angle")
##' calculate_branch_mid(d, layout = "equal_angle")
##' @export
calculate_branch_mid <- function(res, layout) {
    if (layout %in% c("equal_angle", "daylight", "ape", "tree_and_leaf")){
        res$branch.y <- with(res, (y[match(parent, node)] + y)/2)
        res$branch.y[is.na(res$branch.y)] <- 0
    }
    res$branch <- with(res, (x[match(parent, node)] + x)/2)
    if (!is.null(res[['branch.length']])) {
        res$branch.length[is.na(res$branch.length)] <- 0
    }
    res$branch[is.na(res$branch)] <- 0
    if (layout %in% c("equal_angle", "daylight", "ape", "tree_and_leaf")){
        res$branch.x <- res$branch
    }
    return(res)
}


##' Re-assign y coordinates
##'
##' Fills in missing y values bottom-up: every node whose `y` is `NA` gets the
##' mean of its children's y values. Rows with a y value already present are
##' left untouched.
##'
##' @param df a node table with `node`, `parent` and `y` columns.
##' @param currentNode node ids to start the fill from, usually the tips.
##' @return the node table with the missing y values replaced.
##' @examples
##' d <- tree_layout(ape::rtree(6))
##' d$y[c(2, 3)] <- NA
##' re_assign_ycoord_df(d, currentNode = d$node[d$isTip])
##' @export
re_assign_ycoord_df <- function(df, currentNode) {
    if (!anyNA(df$y)) {
        return(df)
    }

    ## One upward pass over the nodes that are reachable from `currentNode`.
    ##
    ## The previous loop advanced one level per iteration and, inside, ran
    ## match()/`==`/`%in%` against the whole table several times. On a
    ## ladder-like tree (depth ~ n) that made this O(n^2): 0.70 s at 300 tips,
    ## 21 s at 1200.
    node <- df$node
    parent <- df$parent
    y <- df$y

    max_id <- max(node, parent, na.rm = TRUE)
    row_of <- rep.int(NA_integer_, max_id)
    row_of[node] <- seq_along(node)

    ## children of every node, by row
    is_child <- !is.na(parent) & parent != node
    kids <- vector("list", nrow(df))
    cl <- split(which(is_child), row_of[parent[is_child]])
    kids[as.integer(names(cl))] <- cl
    kids[vapply(kids, is.null, logical(1))] <- list(integer())

    ## depth within the table, so children are finished before their parent
    depth <- rep.int(NA_integer_, nrow(df))
    roots <- which(!is_child)
    depth[roots] <- 0L
    frontier <- roots
    d <- 0L
    while (length(frontier) > 0L) {
        d <- d + 1L
        nxt <- vector("list", length(frontier))
        for (i in seq_along(frontier)) {
            k <- kids[[frontier[i]]]
            if (length(k) > 0L) {
                depth[k] <- d
                nxt[[i]] <- k
            }
        }
        frontier <- unlist(nxt, use.names = FALSE)
    }

    reachable <- !is.na(match(node, currentNode))
    for (r in order(depth, decreasing = TRUE, na.last = NA)) {
        if (is.na(y[r])) {
            k <- kids[[r]]
            k <- k[!is.na(y[k])]
            if (length(k) > 0L) {
                y[r] <- mean(y[k])
            }
        }
    }

    df$y <- y
    return(df)
}


layoutApe <- function(model, branch.length="branch.length") {
	tree <- as.phylo(model) |> stats::reorder("postorder")

	if (! is.null(tree$edge.length)) {
		if (anyNA(tree$edge.length)) {
			warning("'edge.length' contains NA values...\n## setting 'edge.length' to NULL automatically when plotting the tree...")
			tree$edge.length <- NULL
		}
	}

	if (is.null(tree$edge.length) || branch.length == "none") {
		tree <- set_branch_length_cladogram(tree)
	}

	edge <- tree$edge
	edge.length <- tree$edge.length
	nb.sp <- ape::node.depth(tree)

	df <- as_tibble(model) |>
		mutate(isTip = ! .data$node %in% .data$parent)
	#df$branch.length <- edge.length[df$node] # for cladogram

	# unrooted layout from cran/ape
	M <- ape::unrooted.xy(Ntip(tree),
						  Nnode(tree),
						  tree$edge,
						  tree$edge.length,
						  nb.sp,
						  0)$M
	xx <- M[, 1]
	yy <- M[, 2]

	M <- tibble::tibble(
		node = 1:(Ntip(tree) + Nnode(tree)),
		x = xx - min(xx),
		y = yy - min(yy)
	)

	tree_df <- dplyr::full_join(df, M, by = "node") |>
		as_tibble()
	class(tree_df) <- c("tbl_tree", class(tree_df))
	tree_df
}

.nodeId <- function (tree, type = "all"){
    type <- match.arg(type, c("all", "tips", "internal"))
    if (inherits(tree, "treedata")) {
        tree <- tree@phylo
    }
    nodes <- unique(as.vector(tree$edge))
    if (type == "all") {
        return(nodes)
    }
    edge <- tree$edge
    tips <- edge[!edge[, 2] %in% edge[, 1], 2]
    if (type == "tips"){
        return(tips)
    }
    else if (type == "internal") {
        return(setdiff(nodes, tips))
    }
}







edge2vec <- function(tr) {
  parent <- tr$edge[,1]
  child <- tr$edge[,2]
  
  ## use lookup table
  pvec <- integer(max(tr$edge))
  pvec[child] <- parent
  return(pvec)
}



.as.phylo.hclust2 <- function(x, hang=0.1, ...) {
  h <- x
  tr <- ape::as.phylo(x)
  ev <- edge2vec(tr)

  #extract_inode_hclust_item <- function(h, i, ev) {
  #  j <- h$merge[i,]
  #  if (any(j < 0)) {
  #    j2 <- j[j < 0][1]
  #    res <- ev[abs(j2)]
  #  } else {
  #    res <- ev[extract_inode_hclust_item(h, j, ev)]
  #  }
  #  return(res)
  #}

  #nodes <- vapply(seq_along(h$height), function(i) {
  #  extract_inode_hclust_item(h, i, ev)
  #}, numeric(1))
  
  nodes <- integer(length(h$height))
  for (i in seq_along(nodes)) {
    j <- h$merge[i,]
    if (any(j < 0)) {
      j2 <- j[j < 0][1]
      nodes[i] <- ev[abs(j2)]
    } else {
      nodes[i] <- ev[nodes[j[1]]]
    }
  }

  #len <- numeric(max(tr$edge))
  #len[nodes] <- h$height
  #pn <- ev[nodes]
  #pn[pn == 0] <- tidytree::rootnode(tr)
  #len[nodes] <- len[pn] - len[nodes]
  #len[1:Ntip(tr)] <- hang #max(h$height)/10

  #tr$edge.length <- len[tr$edge[,2]]

  tip2parent <- tr$edge[match(seq_len(Ntip(tr)), tr$edge[,2]), 1]
  if (hang > 0){
    tip.edge.len <- hang * max(h$height) - h$height[match(tip2parent, nodes)]
    attr(tr,'tip.edge.len') <- tip.edge.len
  }
  tr$edge.length <- tr$edge.length * 2
  return(tr)
}