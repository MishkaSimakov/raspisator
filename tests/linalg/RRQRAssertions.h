#pragma once

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "linalg/Linalg.h"

// Shared helpers for verifying the documented rrqr contract (see RRQR.h),
// used by both the dense (RRQRTests.cpp) and sparse (SparseRRQRTests.cpp)
// test suites.
//
//     P A = R Q,   R = [ R1  0 ]   R1 : r x r lower triangular, |R1_ii| > tol,
//                      [ R2 R3 ]   R3 rows have 2-norm <= tol,
//
// where Q is orthogonal and NOT returned. Because Q is orthogonal it cancels
// in the Gram product, giving a Q-free correctness check:
//
//     (P A)(P A)^T = R Q Q^T R^T = R R^T.
//
// That identity, together with the structural checks below, pins the function
// down without ever materialising Q.

namespace rrqr_test {

using Field = double;

inline constexpr Field kEps = 1e-9;

// Row-permuted matrix P A: row i of the result is row perm[i] of A.
inline Matrix<Field> permute_rows(const Matrix<Field>& a,
                                  const std::vector<size_t>& perm) {
  const auto [n, d] = a.shape();
  linalg::Matrix<Field> result(n, d);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < d; ++j) {
      result[i, j] = a[perm[i], j];
    }
  }
  return result;
}

// Gram matrix M M^T.
inline linalg::Matrix<Field> gram(const linalg::Matrix<Field>& m) {
  const auto [n, d] = m.shape();
  linalg::Matrix<Field> g(n, n);
  for (size_t i = 0; i < n; ++i) {
    for (size_t k = 0; k < n; ++k) {
      Field sum = 0;
      for (size_t j = 0; j < d; ++j) {
        sum += m[i, j] * m[k, j];
      }
      g[i, k] = sum;
    }
  }
  return g;
}

inline void expect_matrices_near(const Matrix<Field>& a, const Matrix<Field>& b,
                                 Field eps = kEps) {
  ASSERT_EQ(a.shape(), b.shape());
  const auto [n, d] = a.shape();
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < d; ++j) {
      EXPECT_NEAR((a[i, j]), (b[i, j]), eps) << "at (" << i << ", " << j << ")";
    }
  }
}

inline bool is_permutation(const std::vector<size_t>& perm, size_t n) {
  if (perm.size() != n) {
    return false;
  }
  std::vector<size_t> sorted = perm;
  std::ranges::sort(sorted);
  for (size_t i = 0; i < n; ++i) {
    if (sorted[i] != i) {
      return false;
    }
  }
  return true;
}

// 2-norm of row `row` restricted to columns [from_col, cols).
inline Field row_residual_norm(const linalg::Matrix<Field>& m, size_t row,
                               size_t from_col) {
  Field sum = 0;
  for (size_t j = from_col; j < m.cols(); ++j) {
    sum += m[row, j] * m[row, j];
  }
  return std::sqrt(sum);
}

// Verifies the whole documented rrqr contract, given the original matrix `a`,
// the returned R factor, the row permutation, the numerical rank and the
// tolerance. Independent of how R/perm/rank were obtained, so it serves both
// the dense and the sparse decomposition variants.
inline void verify_rrqr_contract(const linalg::Matrix<Field>& a,
                                 const linalg::Matrix<Field>& r_matrix,
                                 const std::vector<size_t>& perm, size_t rank,
                                 Field tol) {
  const auto [n, d] = a.shape();

  // perm is a genuine permutation of [0, n) ...
  EXPECT_TRUE(is_permutation(perm, n));
  // ... and 0 <= r <= min(n, d).
  EXPECT_LE(rank, std::min(n, d));

  // Q-free correctness: R R^T == (P A)(P A)^T.
  expect_matrices_near(gram(r_matrix), gram(permute_rows(a, perm)));

  // First r rows are lower triangular: everything strictly right of the
  // diagonal vanishes. This covers both "R1 lower triangular" and the zero
  // block in the top-right.
  for (size_t i = 0; i < rank; ++i) {
    for (size_t j = i + 1; j < d; ++j) {
      EXPECT_NEAR((r_matrix[i, j]), 0, kEps)
          << "R1 not lower triangular at (" << i << ", " << j << ")";
    }
  }

  // R1 diagonal: above tolerance and non-increasing in magnitude (the
  // Businger-Golub pivot ordering).
  for (size_t i = 0; i < rank; ++i) {
    EXPECT_GT(std::abs(r_matrix[i, i]), tol) << "tiny pivot at " << i;
  }
  for (size_t i = 0; i + 1 < rank; ++i) {
    EXPECT_GE(std::abs(r_matrix[i, i]) + kEps, std::abs(r_matrix[i + 1, i + 1]))
        << "pivots not non-increasing at " << i;
  }

  // Trailing rows (R3) are linearly dependent: residual 2-norm <= tolerance.
  for (size_t i = rank; i < n; ++i) {
    EXPECT_LE(row_residual_norm(r_matrix, i, rank), tol + kEps)
        << "dependent row " << i << " has residual above tolerance";
  }
}

}  // namespace rrqr_test
