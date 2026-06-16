#pragma once

#include <numeric>
#include <vector>

#include "Matrix.h"
#include "Vector.h"
#include "utils/Accumulators.h"

namespace linalg {

// Rank-revealing QR: column-pivoted Householder QR applied to A^T, used to
// detect a maximal well-conditioned set of linearly independent ROWS of A.
//
// Input : matrix = n x d matrix A (overwritten). `tolerance` is an absolute
//         bound on a row's residual 2-norm.
// Output: matrix is overwritten by R with
//
//             P A = R Q
//
//         P : n x n row permutation,  Q : d x d orthogonal (NOT formed),
//         R : n x d, block form
//
//             [ R1  0 ]   R1 : r x r lower triangular, nonsingular
//             [ R2 R3 ]   0  : r x (d-r) exact zero
//                         R2 : (n-r) x r
//                         R3 : (n-r) x (d-r), each row has 2-norm <= tolerance
//
//         r = numerical rank of A at `tolerance`, 0 <= r <= min(n,d), and
//         |R1_11| >= |R1_22| >= ... >= |R1_rr| > tolerance.
// Return: { perm (size n), r }.  perm[i] = j : row i of P A is row j of A.
//         perm[0..r-1] kept (independent); perm[r..n-1] dependent (removable).
// Thanks to this guide:
// https://hua-zhou.github.io/teaching/biostatm280-2019spring/slides/11-qr/qr.html
template <typename Field>
std::pair<std::vector<size_t>, size_t> rrqr(Matrix<Field>& matrix,
                                            Field tolerance) {
  using std::abs, std::sqrt;

  const auto [n, d] = matrix.shape();
  std::vector<size_t> row_permutation(n);
  std::iota(row_permutation.begin(), row_permutation.end(), 0);

  std::vector<Field> norms_sqr(n, 0);

  for (size_t i = 0; i < n; ++i) {
    Field norm_sqr = 0;

    for (size_t j = 0; j < d; ++j) {
      norms_sqr += matrix[i, j] * matrix[i, j];
    }

    norms_sqr[i] = norm_sqr;
  }

  for (size_t row = 0; row < n; ++row) {
    // choose row with the largest l_2 norm
    ArgMaximum<Field> max_l2;

    for (size_t i = row; i < n; ++i) {
      max_l2.record(i, norms_sqr[i]);
    }

    // stopping criterion
    if (max_l2->max <= tolerance * tolerance) {
      return {std::move(row_permutation), row};
    }

    // swap rows so that pivot row is the current one
    std::swap(row_permutation[row], row_permutation[max_l2->index]);
    std::swap(norms_sqr[row], norms_sqr[max_l2->index]);

    for (size_t i = 0; i < d; ++i) {
      std::swap(matrix[row, i], matrix[max_l2->index, i]);
    }

    const Field norm_sqr = max_l2->max;
    const Field norm = sqrt(norm_sqr);

    // squared norm of u vector from which the Householder matrix is formed:
    // H = I +- 2 * u * u.T / |u|^2
    const Field sign = matrix[row, row] > 0 ? Field(1) : -Field(1);
    const Field new_pivot = matrix[row, row] + sign * norm;
    const Field new_norm_sqr =
        norm_sqr - matrix[row, row] * matrix[row, row] + new_pivot * new_pivot;

    const Field coef_multiplier = Field(2) / new_norm_sqr;

    // go through all rows below and calculate r_i * H_i,
    // where H_i is the Householder matrix for the current row
    for (size_t other_row = row + 1; other_row < n; ++other_row) {
      Field coef = 0;
      for (size_t j = row + 1; j < d; ++j) {
        coef += matrix[row, j] * matrix[other_row, j];
      }
      coef += matrix[other_row, row] * new_pivot;
      coef *= coef_multiplier;

      for (size_t j = row + 1; j < d; ++j) {
        matrix[other_row, j] -= coef * matrix[row, j];
      }
      matrix[other_row, row] -= coef * new_pivot;

      norms_sqr[other_row] -= matrix[other_row, row] * matrix[other_row, row];
    }

    matrix[row, row] = -sign * norm;
    for (size_t i = row + 1; i < d; ++i) {
      matrix[row, i] = 0;
    }
  }

  return {std::move(row_permutation), n};
}

}  // namespace linalg
