// C++ port of the inner loop of .daylightGetArcAngles().
//
// This is a line-by-line transcription, not a rewrite. The R version
// accumulates `arc` in place, and every comparison reads the value written by
// the previous iteration -- note that `arc_left` is compared against
// `theta_child_adj` but stores `theta_child`, which are different numbers when
// the branch wraps around. Reordering the loop or folding the accumulation
// into a max()/min() changes the output. Keep the order.

#include <Rcpp.h>

// std::atan2, std::fabs, std::isnan and std::acos come from <cmath>; Rcpp.h
// happens to pull it in on gcc, but MSVC is not obliged to.
#include <cmath>
#include <vector>

using namespace Rcpp;

// M_PI is not defined by <cmath> on MSVC without _USE_MATH_DEFINES, so compute
// it once instead of relying on the macro.
static const double kPi = std::acos(-1.0);

// [[Rcpp::export]]
NumericVector daylight_arc_angles_cpp(const NumericVector& df_x,
                                      const NumericVector& df_y,
                                      const IntegerVector& row_index,
                                      int origin_id,
                                      const IntegerVector& subtree_node_ids,
                                      const List& children,
                                      const LogicalVector& has_children,
                                      double theta_left_start,
                                      double theta_right_start)
{
    const int origin_row = row_index[origin_id - 1] - 1;
    const double x_origin = df_x[origin_row];
    const double y_origin = df_y[origin_row];

    double theta_left = theta_left_start;
    double theta_right = theta_right_start;

    if (std::isnan(theta_left) || subtree_node_ids.size() == 0) {
        return NumericVector::create(Named("left", 0.0), Named("right", 0.0));
    }

    double arc_left = theta_left;
    double arc_right = theta_right;

    // Angles of every node the loop below reaches, computed once. atan2() is a
    // pure function, so hoisting it is bit-identical to calling it per node.
    const R_xlen_t n = subtree_node_ids.size();
    const std::size_t n_slots = static_cast<std::size_t>(df_x.size());
    std::vector<double> theta_by_id(n_slots, NA_REAL);

    {
        std::vector<int> looked_at;
        looked_at.reserve(static_cast<std::size_t>(n) * 2);
        for (R_xlen_t i = 0; i < n; ++i) {
            const int id = subtree_node_ids[i];
            if (id == origin_id) continue;
            if (!has_children[id - 1]) continue;
            looked_at.push_back(id);
            const IntegerVector ch = children[id - 1];
            for (R_xlen_t j = 0; j < ch.size(); ++j) {
                looked_at.push_back(ch[j]);
            }
        }
        for (std::size_t k = 0; k < looked_at.size(); ++k) {
            const int id = looked_at[k];
            const int row = row_index[id - 1] - 1;
            theta_by_id[static_cast<std::size_t>(id - 1)] =
                atan2(df_y[row] - y_origin, df_x[row] - x_origin) / kPi;
        }
    }

    for (R_xlen_t i = 0; i < n; ++i) {
        const int parent_id = subtree_node_ids[i];
        if (parent_id == origin_id) continue;
        if (!has_children[parent_id - 1]) continue;

        const double theta_parent = theta_by_id[static_cast<std::size_t>(parent_id - 1)];
        const IntegerVector child_ids = children[parent_id - 1];
        if (child_ids.size() == 0) continue;

        for (R_xlen_t j = 0; j < child_ids.size(); ++j) {
            const int child_id = child_ids[j];
            if (child_id == origin_id) continue;

            const double theta_child = theta_by_id[static_cast<std::size_t>(child_id - 1)];

            if ((arc_left < arc_right && !(theta_child > arc_left && theta_child < arc_right)) ||
                (arc_left > arc_right && (theta_child < arc_left && theta_child > arc_right))) {
                continue;
            }

            const double delta = theta_child - theta_parent;
            double delta_adj = delta;
            if (delta > 1) {
                delta_adj = delta - 2;
            } else if (delta < -1) {
                delta_adj = delta + 2;
            }

            double theta_child_adj = theta_child;
            if (delta_adj > 0) {
                if (std::fabs(delta) > 1) {
                    if (arc_left > 0 && theta_child < 0) {
                        theta_child_adj = theta_child + 2;
                    } else if (arc_left < 0 && theta_child > 0) {
                        theta_child_adj = theta_child - 2;
                    }
                }
                // stores theta_child, compares theta_child_adj: intentional
                if (arc_left < theta_child_adj) {
                    arc_left = theta_child;
                }
            } else if (delta_adj < 0) {
                if (std::fabs(delta) > 1) {
                    if (arc_right > 0 && theta_child < 0) {
                        theta_child_adj = theta_child + 2;
                    } else if (arc_right < 0 && theta_child > 0) {
                        theta_child_adj = theta_child - 2;
                    }
                }
                if (arc_right > theta_child_adj) {
                    arc_right = theta_child;
                }
            }
        }
    }

    if (arc_left < 0) arc_left += 2;
    if (arc_right < 0) arc_right += 2;

    return NumericVector::create(Named("left", arc_left), Named("right", arc_right));
}
