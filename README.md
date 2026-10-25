# treelayout: Tree Layout Algorithms

Compute node coordinates for phylogenetic trees and other tree-like
structures, without pulling in a plotting library.

`treelayout` is the layout engine behind
[ggtree](https://github.com/YuLab-SMU/ggtree): it takes a tree and returns
a table with one row per node and an `x`/`y` position for each. Anything
that can draw lines can then draw the tree -- ggplot2, base graphics,
grid, or a downstream format.

## Installation

```r
install.packages("treelayout")

# development version
remotes::install_github("YuLab-SMU/treelayout")
```

## Usage

```r
library(treelayout)
library(ape)

set.seed(42)
tr <- rtree(10)

tree_layout(tr)
#> # A tibble: 19 x 9
#>    parent  node branch.length label isTip     x     y branch angle
#>     <int> <int>         <dbl> <chr> <lgl> <dbl> <dbl>  <dbl> <dbl>
#>  1     11     1         0.560 t10   TRUE  0.560     1  0.280    36
#>  2     13     2         0.989 t1    TRUE  2.03      6  1.54    216
#>  3     14     3         0.082 t8    TRUE  2.07      7  2.03    252
#> ...
#> # i 16 more rows
```

Fifteen layouts are available:

```r
layout_names()
#>  [1] "rectangular"    "slanted"        "fan"            "circular"
#>  [5] "inward_circular" "radial"        "unrooted"       "equal_angle"
#>  [9] "daylight"       "tree_and_leaf"  "dendrogram"     "ape"
#> [13] "ellipse"        "roundrect"      "tidy"
```

```r
tree_layout(tr, layout = "circular")
tree_layout(tr, layout = "tidy")
```

An external layout function works too. It is called with the tree as an
`igraph` object and must return one `x`/`y` row per node:

```r
tree_layout(tr, layout = igraph::layout_with_kk)
```

The layout used is recorded on the result:

```r
get_layout(tree_layout(tr, layout = "circular"))
#> [1] "circular"
```

### Drawing without ggplot2

The point of the package is that the coordinates are just numbers:

```r
d <- tree_layout(tr)

plot(NA, xlim = range(d$x), ylim = range(d$y), asp = 1,
     axes = FALSE, xlab = "", ylab = "")
segments(d$x[match(d$parent, d$node)], d$y[match(d$parent, d$node)],
         d$x, d$y)
text(d$x[d$isTip], d$y[d$isTip], d$label[d$isTip], pos = 4, cex = .7)
```

## Layout families

Not every layout is interchangeable.

| Layout | `x` encodes | tips on evenly spaced `y` |
|---|---|---|
| `rectangular`, `slanted`, `circular`, `fan`, `dendrogram`, `ellipse`, `roundrect`, `tidy` | branch length | yes |
| `equal_angle`, `daylight`, `ape`, `tree_and_leaf` (unrooted) | arbitrary | **no** |
| any external `igraph` layout | arbitrary | **no** |

Only the first group keeps one tip per row, so only it can be combined
with panels that align to each tip (a heat map beside the tree, say).
You can check before drawing:

```r
tip_rows_even <- function(d) {
    y <- sort(d$y[d$isTip])
    if (length(y) < 2) return(TRUE)
    gaps <- diff(y)
    all(gaps > 0) && max(gaps) - min(gaps) < 1e-6
}

tip_rows_even(tree_layout(tr))
#> [1] TRUE
tip_rows_even(tree_layout(tr, layout = "equal_angle"))
#> [1] FALSE
```

## Related packages

| Package | Role |
|---|---|
| [tidytree](https://github.com/YuLab-SMU/tidytree) | the tree data structure this package operates on |
| [treeio](https://github.com/YuLab-SMU/treeio) | reading and writing tree files |
| [ggtree](https://github.com/YuLab-SMU/ggtree) | drawing trees with ggplot2, using these coordinates |

## Author

Guangchuang YU

School of Basic Medical Sciences, Southern Medical University

<https://yulab-smu.top>

## License

Artistic-2.0
