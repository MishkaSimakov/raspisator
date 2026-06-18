#pragma once

#include <numeric>
#include <vector>

#include "CSCMatrix.h"
#include "Matrix.h"
#include "Vector.h"
#include "utils/Accumulators.h"
#include "utils/MergeSortedRange.h"

namespace linalg {

template <typename Field>
struct RRQRResult {
  Matrix<Field> R;
  size_t rank;
  std::vector<size_t> permutation;
};

// Rank-revealing QR: column-pivoted Householder QR applied to A^T, used to
// detect a maximal well-conditioned set of linearly independent ROWS of A.
//
// Input : matrix = n x d matrix A, taken by value (the local copy is consumed
//         and moved into the result). `tolerance` is an absolute bound on a
//         row's residual 2-norm.
// Output: RRQRResult { R, rank, permutation } describing
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
//         rank = r = numerical rank of A at `tolerance`, 0 <= r <= min(n,d),
//         and |R1_11| >= |R1_22| >= ... >= |R1_rr| > tolerance.
//         permutation (size n): permutation[i] = j means row i of P A is row j
//         of A. permutation[0..r-1] are kept (independent); permutation[r..n-1]
//         are dependent (removable).
// Thanks to this guide:
// https://hua-zhou.github.io/teaching/biostatm280-2019spring/slides/11-qr/qr.html
template <typename Field>
RRQRResult<Field> rrqr(Matrix<Field> matrix, Field tolerance) {
  using std::abs, std::sqrt;

  const auto [n, d] = matrix.shape();
  std::vector<size_t> row_permutation(n);
  std::iota(row_permutation.begin(), row_permutation.end(), 0);

  std::vector<Field> norms_sqr(n, 0);

  for (size_t i = 0; i < n; ++i) {
    Field norm_sqr = 0;

    for (size_t j = 0; j < d; ++j) {
      norm_sqr += matrix[i, j] * matrix[i, j];
    }

    norms_sqr[i] = norm_sqr;
  }

  size_t row = 0;
  for (; row < n; ++row) {
    // choose row with the largest l_2 norm
    ArgMaximum<Field> max_l2;

    for (size_t i = row; i < n; ++i) {
      max_l2.record(i, norms_sqr[i]);
    }

    // stopping criterion
    if (max_l2->max <= tolerance * tolerance) {
      break;
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

  return RRQRResult<Field>{
      .R = std::move(matrix),
      .rank = row,
      .permutation = std::move(row_permutation),
  };
}

template <typename Field>
RRQRResult<Field> rrqr(const CSCMatrix<Field>& matrix, Field tolerance) {
  using std::abs, std::sqrt;

  const auto [n, d] = matrix.shape();
  std::vector<size_t> row_permutation(n);
  std::iota(row_permutation.begin(), row_permutation.end(), 0);

  // turn matrix into CSR
  std::vector<std::vector<std::pair<size_t, Field>>> rows(n);
  for (size_t i = 0; i < d; ++i) {
    const size_t col = d - i - 1;

    for (const auto& [row, value] : matrix.get_column(col)) {
      rows[row].emplace_back(col, value);
    }
  }

  // calculate norms
  std::vector<Field> norms_sqr(n, 0);

  for (size_t row = 0; row < n; ++row) {
    Field norm_sqr = 0;

    for (const auto& [col, value] : rows[row]) {
      norm_sqr += value * value;
    }

    norms_sqr[row] = norm_sqr;
  }

  auto result = Matrix<Field>::zeros(n, d);

  // buffer for storing sparse rows in dense format
  std::vector<Field> dense(d, 0);

  std::vector<std::pair<size_t, Field>> sparse_buffer;
  sparse_buffer.reserve(d);

  size_t row = 0;

  for (; row < n; ++row) {
    // choose row with the largest l_2 norm
    ArgMaximum<Field> max_l2;

    for (size_t i = row; i < n; ++i) {
      max_l2.record(i, norms_sqr[i]);
    }

    // stopping criterion
    if (max_l2->max <= tolerance * tolerance) {
      break;
    }

    // ArgMinimum<size_t> min_count;
    //
    // for (size_t i = row; i < n; ++i) {
    //   if (norms_sqr[i] * 10 > max_l2->max) {
    //     min_count.record(i, rows[i].size());
    //   }
    // }

    const size_t pivot_row = max_l2->index;

    // swap rows so that pivot row is the current one
    std::swap(row_permutation[row], row_permutation[pivot_row]);
    std::swap(rows[row], rows[pivot_row]);
    std::swap(norms_sqr[row], norms_sqr[pivot_row]);

    for (size_t i = 0; i < row; ++i) {
      std::swap(result[row, i], result[pivot_row, i]);
    }

    const Field norm_sqr = norms_sqr[row];
    const Field norm = sqrt(norm_sqr);

    // turn current row into dense format
    for (const auto& [col, value] : rows[row]) {
      dense[col] = value;
    }

    // rows[row] is not empty since its norm is greater than tolerance
    assert(!rows[row].empty());
    if (rows[row].back().first == row) {
      rows[row].pop_back();
    }

    // squared norm of u vector from which the Householder matrix is formed:
    // H = I +- 2 * u * u.T / |u|^2
    const Field sign = dense[row] > 0 ? Field(1) : -Field(1);
    const Field new_pivot = dense[row] + sign * norm;
    const Field new_norm_sqr =
        norm_sqr - dense[row] * dense[row] + new_pivot * new_pivot;

    const Field coef_multiplier = Field(2) / new_norm_sqr;

    // go through all rows below and calculate r_i * H_i,
    // where H_i is the Householder matrix for the current row
    for (size_t other_row = row + 1; other_row < n; ++other_row) {
      Field other_row_pivot = 0;

      if (!rows[other_row].empty() && rows[other_row].back().first == row) {
        other_row_pivot = rows[other_row].back().second;
        rows[other_row].pop_back();
      }

      Field coef = 0;
      for (const auto& [col, value] : rows[other_row]) {
        coef += dense[col] * value;
      }
      coef += other_row_pivot * new_pivot;
      coef *= coef_multiplier;

      const Field other_row_new_pivot = other_row_pivot - coef * new_pivot;

      result[other_row, row] = other_row_new_pivot;
      norms_sqr[other_row] -= other_row_new_pivot * other_row_new_pivot;

      if (coef == 0) {
        continue;
      }

      // sparse buffer is empty but has allocated memory in it
      std::swap(rows[other_row], sparse_buffer);

      auto merged =
          MergeSortedRange<Field, std::greater<>>(sparse_buffer, rows[row]);

      for (const auto& [col, left, right] : merged) {
        rows[other_row].emplace_back(col, left - coef * right);
      }

      sparse_buffer.clear();
    }

    // reset current row dense vector
    for (const auto& [col, value] : rows[row]) {
      dense[col] = 0;
    }

    result[row, row] = -sign * norm;
  }

  return RRQRResult<Field>{
      .R = std::move(result),
      .rank = row,
      .permutation = std::move(row_permutation),
  };
}

}  // namespace linalg
