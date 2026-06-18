#include <gtest/gtest.h>

#include <utility>
#include <vector>

#include "linalg/CSCMatrix.h"
#include "linalg/Matrix.h"
#include "linalg/RRQR.h"
#include "linalg/RRQRAssertions.h"

using namespace linalg;
using namespace rrqr_test;

// Tests for the sparse overload
//
//     RRQRResult<Field> rrqr(const CSCMatrix<Field>&, Field tolerance);
//
// It returns the same factorisation P A = R Q as the dense overload, only the
// input is stored column-compressed and R / rank / permutation come back in an
// RRQRResult struct. The contract is identical, so we reuse verify_rrqr_contract
// from RRQRAssertions.h; these tests merely feed it the sparse result and add
// sparsity-specific cases (explicit-zero handling, scattered fill, etc.).

namespace {

struct Decomposition {
  Matrix<Field> r_matrix;
  std::vector<size_t> perm;
  size_t rank;
};

// Builds a CSCMatrix from the dense matrix `a` (dropping nothing -- the sparse
// input is bit-for-bit the same matrix), runs the sparse rrqr and checks the
// full contract against `a`.
Decomposition decompose_and_check(const Matrix<Field>& a, Field tol) {
  // drop_tolerance = 0: keep every structural non-zero so the CSC matrix
  // represents exactly `a`, including tiny entries the rank decision hinges on.
  CSCMatrix<Field> sparse_a(a, Field(0));

  RRQRResult<Field> result = rrqr(sparse_a, tol);

  verify_rrqr_contract(a, result.R, result.permutation, result.rank, tol);

  return {std::move(result.R), std::move(result.permutation), result.rank};
}

}  // namespace

TEST(SparseRRQRTests, Identity) {
  auto dec = decompose_and_check(Matrix<Field>::identity(3), 1e-6);
  EXPECT_EQ(dec.rank, 3u);
}

TEST(SparseRRQRTests, FullRowRankSquare) {
  Matrix<Field> a = {{4, 1, 2}, {1, 5, 3}, {2, 3, 6}};
  auto dec = decompose_and_check(a, 1e-6);
  EXPECT_EQ(dec.rank, 3u);
}

TEST(SparseRRQRTests, WideFullRowRank) {
  // n < d: full row rank, no dependent rows (R3 has zero rows).
  Matrix<Field> a = {{1, 2, 0, 0}, {0, 0, 3, 4}};
  auto dec = decompose_and_check(a, 1e-6);
  EXPECT_EQ(dec.rank, 2u);
}

TEST(SparseRRQRTests, TallReducesToColumnRank) {
  // n > d: rank is capped by the number of columns.
  Matrix<Field> a = {{1, 0}, {0, 1}, {1, 1}, {2, 3}};
  auto dec = decompose_and_check(a, 1e-6);
  EXPECT_EQ(dec.rank, 2u);
}

TEST(SparseRRQRTests, ExactLinearDependency) {
  // row 2 = row 0 + row 1, so the numerical rank is 2.
  Matrix<Field> a = {{1, 1, 0, 0}, {0, 0, 1, 1}, {1, 1, 1, 1}};
  auto dec = decompose_and_check(a, 1e-6);
  EXPECT_EQ(dec.rank, 2u);
}

TEST(SparseRRQRTests, DuplicateRows) {
  Matrix<Field> a = {{1, 2, 3}, {1, 2, 3}};
  auto dec = decompose_and_check(a, 1e-6);
  ASSERT_EQ(dec.rank, 1u);
  // One of the two identical rows is kept, the other dropped.
  EXPECT_NE(dec.perm[0], dec.perm[1]);
}

TEST(SparseRRQRTests, ZeroRowIsDropped) {
  Matrix<Field> a = {{1, 0}, {0, 0}, {0, 1}};
  auto dec = decompose_and_check(a, 1e-6);
  ASSERT_EQ(dec.rank, 2u);
  // The zero row (original index 1) is linearly dependent -> ends up last.
  EXPECT_EQ(dec.perm[2], 1u);
}

TEST(SparseRRQRTests, ZeroMatrix) {
  auto dec = decompose_and_check(Matrix<Field>::zeros(3, 3), 1e-6);
  EXPECT_EQ(dec.rank, 0u);
  // No pivot step succeeds, so the permutation is untouched.
  EXPECT_EQ(dec.perm, (std::vector<size_t>{0, 1, 2}));
}

TEST(SparseRRQRTests, EmptyColumns) {
  // A column that is entirely zero contributes no entries to any CSC column;
  // the rank is still governed by the populated columns.
  Matrix<Field> a = {{1, 0, 2}, {0, 0, 0}, {3, 0, 4}};
  auto dec = decompose_and_check(a, 1e-6);
  EXPECT_EQ(dec.rank, 2u);
  // The all-zero middle row is dependent and sorted to the back.
  EXPECT_EQ(dec.perm[2], 1u);
}

TEST(SparseRRQRTests, NumericalRankDependsOnTolerance) {
  // The two rows differ only by a tiny epsilon component, so the second row's
  // residual after orthogonalisation is ~= epsilon. Whether it counts as
  // dependent is decided purely by the tolerance. Built with drop_tolerance 0
  // so the epsilon entry survives the dense -> sparse conversion.
  const Field epsilon = 1e-4;
  Matrix<Field> a = {{1, 0}, {1, epsilon}};

  // tol < epsilon: the rows are numerically independent.
  EXPECT_EQ(decompose_and_check(a, 1e-6).rank, 2u);

  // tol > epsilon: the second row is numerically dependent.
  EXPECT_EQ(decompose_and_check(a, 1e-2).rank, 1u);
}

TEST(SparseRRQRTests, DiagonalIsFullRank) {
  // A genuinely sparse, well-conditioned matrix: every row is independent.
  Matrix<Field> a = {
      {5, 0, 0, 0, 0},
      {0, 4, 0, 0, 0},
      {0, 0, 3, 0, 0},
      {0, 0, 0, 2, 0},
      {0, 0, 0, 0, 1},
  };
  auto dec = decompose_and_check(a, 1e-6);
  EXPECT_EQ(dec.rank, 5u);
}

TEST(SparseRRQRTests, ScatteredSparseWithDependency) {
  // Larger, scattered sparsity pattern with one dependent row:
  // row 4 = row 0 + row 2, so the numerical rank is 4.
  Matrix<Field> a = {
      {1, 0, 0, 2, 0, 0},
      {0, 3, 0, 0, 0, 1},
      {0, 0, 4, 0, 5, 0},
      {0, 0, 0, 0, 0, 7},
      {1, 0, 4, 2, 5, 0},
  };
  auto dec = decompose_and_check(a, 1e-6);
  EXPECT_EQ(dec.rank, 4u);
  // The dependent row (original index 4) is removable -> last in perm.
  EXPECT_EQ(dec.perm.back(), 4u);
}

TEST(SparseRRQRTests, AllRowsDependentOnOne) {
  // Three rows that are scalar multiples of a single sparse direction: rank 1.
  Matrix<Field> a = {{0, 2, 0, 0}, {0, 4, 0, 0}, {0, -1, 0, 0}};
  auto dec = decompose_and_check(a, 1e-6);
  EXPECT_EQ(dec.rank, 1u);
}

TEST(SparseRRQRTests, SingleRow) {
  Matrix<Field> a = {{0, 0, 7, 0}};
  auto dec = decompose_and_check(a, 1e-6);
  EXPECT_EQ(dec.rank, 1u);
}

TEST(SparseRRQRTests, SingleColumn) {
  // n > d = 1: rank capped at 1, the remaining rows are dependent.
  Matrix<Field> a = {{2}, {0}, {3}};
  auto dec = decompose_and_check(a, 1e-6);
  EXPECT_EQ(dec.rank, 1u);
}
