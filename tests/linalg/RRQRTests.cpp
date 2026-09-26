#include <gtest/gtest.h>

#include <algorithm>
#include <utility>
#include <vector>

#include "linalg/Matrix.h"
#include "linalg/RRQR.h"
#include "linalg/RRQRAssertions.h"

using namespace linalg;
using namespace rrqr_test;

// These tests encode the documented contract of rrqr (see RRQR.h). The shared
// verification logic lives in RRQRAssertions.h, which is reused by the sparse
// suite (SparseRRQRTests.cpp).

namespace {

struct Decomposition {
  Matrix<Field> r_matrix;
  std::vector<size_t> perm;
  size_t rank;
};

// Runs the dense rrqr on `a` (passed by value, so the original is left intact)
// and verifies the whole contract against it. Returns the result so individual
// tests can assert the rank and permutation specifics on top.
Decomposition decompose_and_check(const Matrix<Field>& a, Field tol) {
  RRQRResult<Field> result = rrqr(a, tol);

  verify_rrqr_contract(a, result.R, result.permutation, result.rank, tol);

  return {std::move(result.R), std::move(result.permutation), result.rank};
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
