#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#include "linalg/Matrix.h"
#include "linalg/RRQR.h"

using namespace linalg;
using Field = double;

// These tests encode the documented contract of rrqr (see RRQR.h):
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

namespace {

constexpr Field kEps = 1e-9;

// Row-permuted matrix P A: row i of the result is row perm[i] of A.
Matrix<Field> permute_rows(const Matrix<Field>& a,
                           const std::vector<size_t>& perm) {
  const auto [n, d] = a.shape();
  Matrix<Field> result(n, d);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < d; ++j) {
      result[i, j] = a[perm[i], j];
    }
  }
  return result;
}

// Gram matrix M M^T.
Matrix<Field> gram(const Matrix<Field>& m) {
  const auto [n, d] = m.shape();
  Matrix<Field> g(n, n);
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

void expect_matrices_near(const Matrix<Field>& a, const Matrix<Field>& b,
                          Field eps = kEps) {
  ASSERT_EQ(a.shape(), b.shape());
  const auto [n, d] = a.shape();
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < d; ++j) {
      EXPECT_NEAR((a[i, j]), (b[i, j]), eps) << "at (" << i << ", " << j << ")";
    }
  }
}

bool is_permutation(const std::vector<size_t>& perm, size_t n) {
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
Field row_residual_norm(const Matrix<Field>& m, size_t row, size_t from_col) {
  Field sum = 0;
  for (size_t j = from_col; j < m.cols(); ++j) {
    sum += m[row, j] * m[row, j];
  }
  return std::sqrt(sum);
}

struct Decomposition {
  Matrix<Field> r_matrix;
  std::vector<size_t> perm;
  size_t rank;
};

// Runs rrqr on a copy of `a` and verifies the whole contract against the
// original. Returns the result so individual tests can assert the rank and
// permutation specifics on top.
Decomposition decompose_and_check(const Matrix<Field>& a, Field tol) {
  const auto [n, d] = a.shape();
  Matrix<Field> r_matrix = a;
  auto [perm, rank] = rrqr(r_matrix, tol);

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

  return {std::move(r_matrix), std::move(perm), rank};
}

}  // namespace

TEST(RRQRTests, Identity) {
  auto dec = decompose_and_check(Matrix<Field>::identity(3), 1e-6);
  EXPECT_EQ(dec.rank, 3u);
}

TEST(RRQRTests, FullRowRankSquare) {
  Matrix<Field> a = {{4, 1, 2}, {1, 5, 3}, {2, 3, 6}};
  auto dec = decompose_and_check(a, 1e-6);
  EXPECT_EQ(dec.rank, 3u);
}

TEST(RRQRTests, WideFullRowRank) {
  // n < d: full row rank, no dependent rows (R3 has zero rows).
  Matrix<Field> a = {{1, 2, 0, 0}, {0, 0, 3, 4}};
  auto dec = decompose_and_check(a, 1e-6);
  EXPECT_EQ(dec.rank, 2u);
}

TEST(RRQRTests, TallReducesToColumnRank) {
  // n > d: rank is capped by the number of columns.
  Matrix<Field> a = {{1, 0}, {0, 1}, {1, 1}, {2, 3}};
  auto dec = decompose_and_check(a, 1e-6);
  EXPECT_EQ(dec.rank, 2u);
}

TEST(RRQRTests, ExactLinearDependency) {
  // row 2 = row 0 + row 1, so the numerical rank is 2.
  Matrix<Field> a = {{1, 1, 0, 0}, {0, 0, 1, 1}, {1, 1, 1, 1}};
  auto dec = decompose_and_check(a, 1e-6);
  EXPECT_EQ(dec.rank, 2u);
}

TEST(RRQRTests, DuplicateRows) {
  Matrix<Field> a = {{1, 2, 3}, {1, 2, 3}};
  auto dec = decompose_and_check(a, 1e-6);
  ASSERT_EQ(dec.rank, 1u);
  // One of the two identical rows is kept, the other dropped.
  EXPECT_NE(dec.perm[0], dec.perm[1]);
}

TEST(RRQRTests, ZeroRowIsDropped) {
  Matrix<Field> a = {{1, 0}, {0, 0}, {0, 1}};
  auto dec = decompose_and_check(a, 1e-6);
  ASSERT_EQ(dec.rank, 2u);
  // The zero row (original index 1) is linearly dependent -> ends up last.
  EXPECT_EQ(dec.perm[2], 1u);
}

TEST(RRQRTests, ZeroMatrix) {
  auto dec = decompose_and_check(Matrix<Field>::zeros(3, 3), 1e-6);
  EXPECT_EQ(dec.rank, 0u);
  // No pivot step succeeds, so the permutation is untouched.
  EXPECT_EQ(dec.perm, (std::vector<size_t>{0, 1, 2}));
}

TEST(RRQRTests, NumericalRankDependsOnTolerance) {
  // The two rows differ only by a tiny epsilon component, so the second row's
  // residual after orthogonalisation is ~= epsilon. Whether it counts as
  // dependent is decided purely by the tolerance.
  const Field epsilon = 1e-4;
  Matrix<Field> a = {{1, 0}, {1, epsilon}};

  // tol < epsilon: the rows are numerically independent.
  EXPECT_EQ(decompose_and_check(a, 1e-6).rank, 2u);

  // tol > epsilon: the second row is numerically dependent.
  EXPECT_EQ(decompose_and_check(a, 1e-2).rank, 1u);
}
