#' @keywords internal
"_PACKAGE"

## tidytree supplies the traversal generics and `as_tibble()` dispatch that
## turn a tree into a node table. These are imported (not called via `::`)
## so the S3 methods registered inside tidytree can reach their own helpers;
## `tidytree::offspring()` fails at runtime otherwise.
##
## Only what is used is imported. This package never reads treedata slots and
## deliberately does not import the class.
#' @importFrom tidytree child parent ancestor rootnode isTip offspring sibling
#' @importFrom tidytree as_tibble as.phylo Nnode Ntip getNodeNum
#' @importFrom rlang .data
#' @importFrom methods is
NULL
