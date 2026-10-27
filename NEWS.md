# treelayout 0.1.0.003

+ `calculate_angle()` stays within one turn after `scaleClade()` shifts `y`; it used to pass 360 and send tip labels to the wrong side (2026-10-27, #647)
+ remove quadratic paths that made deep trees slow in `getYcoord()`, `getXcoord_no_length()`, `getXcoord2()`, `re_assign_ycoord_df()`, the `yscale` helpers and `slanted`: one post-order pass instead of one level per iteration, 1600 ladder tips from 17.7 s to 0.035 s (2026-10-27)
+ `layout = "unrooted"` runs the daylight method; it fell through to the rectangular branch and returned its coordinates while reporting `"unrooted"` (2026-10-27, #568)
+ `daylight` and `tree_and_leaf` have compiled inner loops, same coordinates as the R reference, 5.8x and 2.9x at 1000 tips; the R versions remain under `options(treelayout.daylight_cpp = FALSE)` (2026-10-27)
+ `tree_layout()` sets `tbl_tree` once in the class vector, not twice (2026-10-27)
+ a layout function returning the wrong shape is reported by name instead of failing with `subscript out of bounds` (2026-10-27)
+ the layouts are iterative; a 3000-tip ladder tree is now a test, guarding the C-stack failure of ggtree#400 (2026-10-27)
+ first CRAN release: the layout engine extracted from [ggtree](https://github.com/YuLab-SMU/ggtree), so tree coordinates can be computed without a plotting dependency (2026-10-27)

## Layouts

Fifteen: `rectangular`, `slanted`, `fan`, `circular`, `inward_circular`,
`radial`, `unrooted`, `equal_angle`, `daylight`, `tree_and_leaf`, `dendrogram`,
`ape`, `ellipse`, `roundrect`, `tidy`. `unrooted` is an alias for `daylight`.

A function can also be passed as `layout`; it receives the tree as an `igraph`
object and must return one `x`/`y` row per node. `igraph` is in `Suggests`.

## Notes

+ Imports `ape`, `tidytree`, `tibble`, `dplyr`, `rlang`, `stats`, `methods`,
  `utils` and `Rcpp`.
+ Coordinates are identical to ggtree 4.3.1 for all fifteen layouts.
+ `tree_and_leaf` makes every tip repel every other on every iteration, so its
  cost is quadratic in the number of tips by construction.
