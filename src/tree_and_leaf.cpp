// C++ port of the force-directed relaxation in layoutTreeAndLeaf().
//
// Line-by-line transcription of the R loop, in the same order. Unlike the
// daylight arc accumulation, every quantity here is either element-wise or a
// sum of contributions, so the arithmetic is order-independent for identical
// iteration counts; the loop order itself is still preserved exactly.
//
// Fallback directions for coincident points reproduce
// .treeAndLeafFallbackDirections(): angle = ((id_a + id_b) * 0.6180339887...)
// mod 1 * 2 * pi, which must match R's fmod semantics for negative sums.

#include <Rcpp.h>
#include <cmath>
#include <vector>

// Rcpp's Matrix/Vector operator() takes a signed index type; indexing with int
// is the idiomatic form across every Rcpp package. -Wsign-conversion flags the
// implicit int -> size_t inside the STL containers anyway, and silencing it
// here reads better than casting at 40 call sites. --as-cran does not enable
// -Wsign-conversion, so this only affects developers compiling with -Wextra.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"

using namespace Rcpp;

// M_PI is not guaranteed by MSVC; <cmath> may not define it without
// _USE_MATH_DEFINES.
static const double kPi = std::acos(-1.0);

// [[Rcpp::export]]
NumericMatrix tree_and_leaf_relax_cpp(NumericMatrix positions,
                                      const NumericVector& depth_norm,
                                      const IntegerVector& tip_rows,      // 1-based
                                      const IntegerVector& internal_rows, // 1-based
                                      int root_row,                       // 1-based
                                      const IntegerVector& edge_parent_rows, // 1-based
                                      const IntegerVector& edge_child_rows,  // 1-based
                                      const NumericVector& target_edge_length,
                                      const NumericVector& outward_dir,   // n x 2 flattened
                                      const NumericVector& mobility,
                                      const IntegerVector& node_ids,      // 1-based ids by row
                                      double leaf_force,
                                      double edge_force,
                                      double anchor_force,
                                      double internal_anchor,
                                      double outward_force,
                                      double step,
                                      double cooling,
                                      double tol,
                                      int max_iter)
{
    const int n_nodes = positions.nrow();
    NumericMatrix initial_positions = clone(positions);
    double current_step = step;

    // tip weights are constant across iterations
    const int n_tips = static_cast<int>(tip_rows.size());
    std::vector<double> tip_w(static_cast<std::size_t>(n_tips));
    for (int i = 0; i < n_tips; ++i)
        tip_w[static_cast<std::size_t>(i)] = 0.5 + depth_norm[tip_rows[i] - 1];

    // anchor weights are constant across iterations
    std::vector<double> anchor_w(static_cast<std::size_t>(n_nodes), anchor_force);
    for (R_xlen_t i = 0; i < internal_rows.size(); ++i)
        anchor_w[static_cast<std::size_t>(internal_rows[i] - 1)] = internal_anchor;
    anchor_w[static_cast<std::size_t>(root_row - 1)] = 0.0;

    const double golden = 0.618033988749895;

    for (int iter = 0; iter < max_iter; ++iter) {
        NumericMatrix forces(n_nodes, 2);
        std::fill(forces.begin(), forces.end(), 0.0);

        // --- tip-tip repulsion ------------------------------------------------
        for (int i = 0; i < n_tips - 1; ++i) {
            const int row_i = tip_rows[i] - 1;
            for (int j = i + 1; j < n_tips; ++j) {
                const int row_j = tip_rows[j] - 1;
                double dx = positions(row_j, 0) - positions(row_i, 0);
                double dy = positions(row_j, 1) - positions(row_i, 1);
                double dist2 = dx * dx + dy * dy;
                if (dist2 < 1e-12) {
                    // fallback direction from the node ids
                    const double angle =
                        std::fmod(((double)node_ids[row_i] + (double)node_ids[row_j]) * golden,
                                  1.0) * 2.0 * kPi;
                    dx = std::cos(angle);
                    dy = std::sin(angle);
                    dist2 = dx * dx + dy * dy;
                }
                const double dist = std::sqrt(dist2);
                const double magnitude =
                    leaf_force * (tip_w[i] + tip_w[j]) / (dist2 + 1e-6);
                const double cx = dx / dist * magnitude;
                const double cy = dy / dist * magnitude;
                forces(row_i, 0) -= cx;
                forces(row_i, 1) -= cy;
                forces(row_j, 0) += cx;
                forces(row_j, 1) += cy;
            }
        }

        // --- edge springs -----------------------------------------------------
        for (int e = 0; e < edge_child_rows.size(); ++e) {
            const int pr = edge_parent_rows[e] - 1;
            const int cr = edge_child_rows[e] - 1;
            double dx = positions(cr, 0) - positions(pr, 0);
            double dy = positions(cr, 1) - positions(pr, 1);
            double dist2 = dx * dx + dy * dy;
            if (dist2 < 1e-12) {
                const double angle =
                    std::fmod(((double)node_ids[pr] + (double)node_ids[cr]) * golden,
                              1.0) * 2.0 * kPi;
                dx = std::cos(angle);
                dy = std::sin(angle);
                dist2 = dx * dx + dy * dy;
            }
            const double dist = std::sqrt(dist2);
            const double magnitude =
                edge_force * (dist - target_edge_length[e]);
            const double cx = dx / dist * magnitude;
            const double cy = dy / dist * magnitude;
            forces(pr, 0) += cx;
            forces(pr, 1) += cy;
            forces(cr, 0) -= cx;
            forces(cr, 1) -= cy;
        }

        // --- anchors ----------------------------------------------------------
        for (int r = 0; r < n_nodes; ++r) {
            forces(r, 0) += (initial_positions(r, 0) - positions(r, 0)) * anchor_w[r];
            forces(r, 1) += (initial_positions(r, 1) - positions(r, 1)) * anchor_w[r];
        }

        // --- outward force on tips -------------------------------------------
        if (outward_force != 0.0) {
            for (int i = 0; i < n_tips; ++i) {
                const int r = tip_rows[i] - 1;
                const double scale =
                    outward_force * (0.5 + depth_norm[r]);
                forces(r, 0) += outward_dir[r * 2 + 0] * scale;
                forces(r, 1) += outward_dir[r * 2 + 1] * scale;
            }
        }

        // --- displacement, clamped -------------------------------------------
        double max_disp = 0.0;
        for (int r = 0; r < n_nodes; ++r) {
            double dx = forces(r, 0) * (current_step * mobility[r]);
            double dy = forces(r, 1) * (current_step * mobility[r]);
            const double norm = std::sqrt(dx * dx + dy * dy);
            if (norm > current_step && norm > 0.0) {
                const double f = current_step / norm;
                dx *= f;
                dy *= f;
            }
            positions(r, 0) += dx;
            positions(r, 1) += dy;
            if (norm > max_disp) max_disp = norm;
        }
        positions(root_row - 1, 0) = initial_positions(root_row - 1, 0);
        positions(root_row - 1, 1) = initial_positions(root_row - 1, 1);

        if (max_disp < tol) break;
        current_step *= cooling;
    }

    return positions;
}

#pragma GCC diagnostic pop
