# treelayout 0.1.0.004

+ `calculate_angle()` stays within one turn after `scaleClade()` shifts `y`; it used to pass 360 and send tip labels to the wrong side (2026-10-27, ggtree#647)
+ remove quadratic paths that made deep trees slow in `getYcoord()`, `getXcoord_no_length()`, `getXcoord2()`, `re_assign_ycoord_df()`, the `yscale` helpers and `slanted`: one post-order pass instead of one level per iteration, 1600 ladder tips from 17.7 s to 0.035 s (2026-10-27)
+ `layout = "unrooted"` runs the daylight method; it fell through to the rectangular branch and returned its coordinates while reporting `"unrooted"` (2026-10-27, ggtree#568)
+ `daylight` and `tree_and_leaf` have compiled inner loops, same coordinates as the R reference, 5.8x and 2.9x at 1000 tips; the R versions remain under `options(treelayout.daylight_cpp = FALSE)` (2026-10-27)
+ `tree_layout()` sets `tbl_tree` once in the class vector, not twice (2026-10-27)
+ a layout function returning the wrong shape is reported by name instead of failing with `subscript out of bounds` (2026-10-27)
+ the layouts are iterative; a 3000-tip ladder tree is now a test, guarding the C-stack failure of ggtree#400 (2026-10-27)

# treelayout 0.1.0

First CRAN release. `treelayout` is the layout engine extracted from
[ggtree](https://github.com/YuLab-SMU/ggtree), so that tree coordinates can be
computed without a plotting dependency.

## Public API

+ `tree_layout()` turns a tree into a table of node coordinates, one row per
  node. It handles `phylo`, `treedata`, `multiPhylo` and `hclust`-like objects
  (`hclust`, `dendrogram`, `agnes`, `diana`, `twins`, `linkage`).
+ `get_layout()` reads back the layout name recorded on a result.
+ `layout_names()` lists the layouts `tree_layout()` accepts by name.
+ `clade_geometry()` returns the arc angles and radii of a clade, for code that
  annotates clades.
+ `Date2decimal()` / `decimal2Date()` convert between dates and decimal years,
  for time-scaled trees.
+ `getXcoord2()`, `getYcoord()`, `calculate_angle()`,
  `calculate_branch_mid()`, `re_assign_ycoord_df()` are layout primitives for
  code that changes coordinates after the initial layout and has to recompute
  the derived columns.

## Layouts

Fifteen: `rectangular`, `slanted`, `fan`, `circular`, `inward_circular`,
`radial`, `unrooted`, `equal_angle`, `daylight`, `tree_and_leaf`, `dendrogram`,
`ape`, `ellipse`, `roundrect`, `tidy`.

`unrooted` is an alias for `daylight`, the name ggtree users know; `get_layout()`
records the method that actually ran, `"daylight"`.

A function can be passed as `layout` as well. It is called with the tree
converted to an `igraph` object and must return one `x`/`y` row per node, so
any layout from igraph or graphlayouts can be used.

## Compiled hot paths

`daylight` and `tree_and_leaf` are the only layouts that iterate. Their inner
loops are compiled (Rcpp) and give the same coordinates as the R reference.
Measured single-threaded at 1000 tips, each run in a fresh R process:

+ `daylight` 33.0 s -> 5.7 s, about 5.8x
+ `tree_and_leaf` 52.7 s -> 18.4 s, about 2.9x

The R implementations remain under option `treelayout.daylight_cpp = FALSE` and
the tests compare the two on every platform CI covers.

`tree_and_leaf` is a force-directed relaxation: every tip repels every other on
every iteration, so its cost is quadratic in the number of tips by construction.
That is why its speedup is smaller than `daylight`'s -- the compiled loop
removes interpreter overhead, not the pairwise work.

## Quadratic paths removed

Several loops advanced one tree level per iteration, so their cost grew with
tree depth rather than with the number of nodes. On a ladder-like tree, where
depth is about the tip count, that is quadratic. They now walk a single
post-order pass (`order(ape::node.depth(tree))`, which visits every child before
its parent):

+ `getYcoord()` and `getXcoord_no_length()`, reached by every layout that is not
  unrooted, tidy or external: 1600 ladder tips went from 17.7 s to 0.035 s.
+ `getXcoord2()`, used by `slanted`: 2000 ladder tips, 4.74 s -> 0.038 s.
+ `re_assign_ycoord_df()`: 1200 ladder tips, 21.2 s -> 0.010 s.
+ The `yscale` helpers called `child()` / `parent()` once per node, each of which
  scans the whole edge list: 2000 tips, 5.4 s -> 0.067 s.
+ `slanted` built an `n_tip x n_node` indicator matrix and then summarised it.
  The summaries now accumulate in one pass: 2000 ladder tips, 4.8 s -> 0.038 s.

`daylight` also reuses its angle scratch buffer instead of allocating one per
subtree (46 MB per round at 1000 tips), and runs its rotation loop on a plain
data frame rather than a `tbl_df`.

## Notes on the extraction

+ No plotting dependency: the package imports only `ape`, `tidytree`, `tibble`,
  `dplyr`, `rlang`, `stats`, `methods`, `utils` and `Rcpp`. `igraph` is in
  `Suggests` and is only needed when a layout function is passed.
+ Coordinates are identical to ggtree 4.3.1 for all fifteen layouts, verified in
  `ggtree`'s own test suite.
+ `layoutDaylight()` runs up to `MAX_COUNT` rounds and stops when the average
  angle change drops below `MINIMUM_AVERAGE_ANGLE_CHANGE`. An attempt to also
  stop as soon as a round stopped improving was reverted: on trees that keep
  oscillating it moved the geometry by up to the full x range, and those trees
  were not converging anyway.
+ The per-round "Average angle change" message from `layoutDaylight()` is behind
  `getOption("treelayout.verbose")`, so layouts are quiet by default.
+ The internal functions carry roxygen comments but `@noRd`, so they stay
  documented in the source without producing help pages.
