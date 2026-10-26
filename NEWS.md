# treelayout 0.0.0.9000

First release. `treelayout` is the layout engine extracted from
[ggtree](https://github.com/YuLab-SMU/ggtree), so that tree coordinates can be
computed without a plotting dependency.

## Public API

+ `tree_layout()` turns a tree into a table of node coordinates, one row per
  node. It handles `phylo`, `treedata`, `multiPhylo` and `hclust`-like objects.
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

## Compiled hot paths

`daylight` and `tree_and_leaf` are the only layouts that iterate. Their inner
loops are compiled (Rcpp) and give the same coordinates as the R reference:

+ `daylight` 3.9-120s -> 1.2-18s, about 5.7x faster at 1000 tips
+ `tree_and_leaf` 8.2-52s -> 1.2-7.1s, about 7.4x faster at 1000 tips

The R implementations remain under option `treelayout.daylight_cpp = FALSE`
and the tests compare the two on every platform CI covers.

## Layouts

Fifteen: `rectangular`, `slanted`, `fan`, `circular`, `inward_circular`,
`radial`, `unrooted`, `equal_angle`, `daylight`, `tree_and_leaf`, `dendrogram`,
`ape`, `ellipse`, `roundrect`, `tidy`.

A function can be passed as `layout` as well. It is called with the tree
converted to an `igraph` object and must return one `x`/`y` row per node, so
any layout from igraph or graphlayouts can be used.

## Notes on the extraction

+ No plotting dependency: the package imports only `ape`, `tidytree`, `tibble`,
  `dplyr`, `purrr`, `rlang`, `stats`, `methods` and `utils`. `igraph` is in
  `Suggests` and is only needed when a layout function is passed.
+ Coordinates are identical to ggtree 4.3.1 for all fifteen layouts, verified
  in `ggtree`'s own test suite.
+ `layoutDaylight()` runs up to `MAX_COUNT` rounds and stops when the average
  angle change drops below `MINIMUM_AVERAGE_ANGLE_CHANGE`. An attempt to also
  stop as soon as a round stopped improving was reverted: on trees that keep
  oscillating it moved the geometry by up to the full x range, and those trees
  were not converging anyway.
+ The per-round "Average angle change" message from `layoutDaylight()` is now
  behind `getOption("treelayout.verbose")`, so layouts are quiet by default.
