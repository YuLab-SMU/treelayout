## Local aliases for tidytree generics used throughout the layout code.
## `getRoot` and `get.tree` keep the names the moved code already uses.

#' @noRd
getRoot <- tidytree::rootnode

#' @noRd
get.tree <- tidytree::as.phylo

## S3 methods for tbl_tree are not exported by tidytree; the layout code calls
## them directly on a node table, so bind the names here.
#' @noRd
.child.tbl_tree <- utils::getFromNamespace("child.tbl_tree", "tidytree")

#' @noRd
.parent.tbl_tree <- utils::getFromNamespace("parent.tbl_tree", "tidytree")
