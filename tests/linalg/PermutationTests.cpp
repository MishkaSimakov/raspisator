#include <gtest/gtest.h>

#include <algorithm>
#include <format>
#include <numeric>
#include <random>
#include <utility>
#include <vector>

#include "ConstructSparse.h"
#include "field/BigInteger.h"
#include "field/FieldTraits.h"
#include "linalg/Linalg.h"
#include "linalg/Permutation.h"

using namespace linalg;

namespace {

// The reference permutation used by the deterministic tests:
//   P[0] = 2, P[1] = 0, P[2] = 3, P[3] = 1
// (single 4-cycle 0 -> 2 -> 3 -> 1 -> 0, an odd permutation)
Permutation reference() { return Permutation::from_vector({2, 0, 3, 1}); }

// A 4x2 matrix with pairwise distinct, easily traceable rows.
Matrix<Rational> rows_4x2() {
  return {
      {10, 11},
      {20, 21},
      {30, 31},
      {40, 41},
  };
}

}  // namespace

// ---------------------------------------------------------------------------
// id
// ---------------------------------------------------------------------------

TEST(PermutationTests, IdHasGivenSize) {
  const auto p = Permutation::id(5);
  ASSERT_EQ(p.size(), 5u);
}

TEST(PermutationTests, IdMapsEveryIndexToItself) {
  const auto p = Permutation::id(5);
  for (size_t i = 0; i < p.size(); ++i) {
    ASSERT_EQ(p[i], i);
  }
}

TEST(PermutationTests, IdLeavesMatrixUnchanged) {
  const auto A = rows_4x2();
  const auto p = Permutation::id(4);
  ASSERT_EQ(p.apply(A), A);
  ASSERT_EQ(p.apply_transposed(A), A);
}

TEST(PermutationTests, IdIsEven) { ASSERT_TRUE(Permutation::id(5).is_even()); }

// ---------------------------------------------------------------------------
// from_vector / size / operator[]
// ---------------------------------------------------------------------------

TEST(PermutationTests, FromVectorPreservesSize) {
  ASSERT_EQ(reference().size(), 4u);
}

TEST(PermutationTests, OperatorIndexReturnsMapping) {
  const auto p = reference();
  ASSERT_EQ(p[0], 2u);
  ASSERT_EQ(p[1], 0u);
  ASSERT_EQ(p[2], 3u);
  ASSERT_EQ(p[3], 1u);
}

TEST(PermutationTests, FromVectorEmptyPermutation) {
  const auto p = Permutation::from_vector({});
  ASSERT_EQ(p.size(), 0u);
}

// ---------------------------------------------------------------------------
// apply (dense Matrix) : result = P A, i.e. row i of A becomes row P[i]
// ---------------------------------------------------------------------------

TEST(PermutationTests, ApplyDensePermutesRows) {
  const auto p = reference();
  const Matrix<Rational> expected = {
      {20, 21},  // row 0 of PA = A row 1 (P[1] = 0)
      {40, 41},  // row 1 of PA = A row 3 (P[3] = 1)
      {10, 11},  // row 2 of PA = A row 0 (P[0] = 2)
      {30, 31},  // row 3 of PA = A row 2 (P[2] = 3)
  };
  ASSERT_EQ(p.apply(rows_4x2()), expected);
}

// ---------------------------------------------------------------------------
// apply (CSCMatrix) : same row permutation, sparse representation
// ---------------------------------------------------------------------------

TEST(PermutationTests, ApplySparseMatchesDense) {
  const auto p = reference();

  const auto A_dense = rows_4x2();
  const auto A_sparse = CSCMatrix(A_dense);

  const Matrix<Rational> from_sparse(p.apply(A_sparse));
  const Matrix<Rational> from_dense = p.apply(A_dense);

  ASSERT_EQ(from_sparse, from_dense);
}

// ---------------------------------------------------------------------------
// apply (vector of (row, value)) : relabels the row index of each entry
// ---------------------------------------------------------------------------

TEST(PermutationTests, ApplySparseVectorRelabelsRows) {
  const auto p = reference();

  std::vector<std::pair<size_t, Rational>> column = {
      {0, Rational{5}},
      {1, Rational{6}},
      {2, Rational{7}},
      {3, Rational{8}},
  };

  const auto result = p.apply(column);

  const std::vector<std::pair<size_t, Rational>> expected = {
      {2, Rational{5}},  // P[0] = 2
      {0, Rational{6}},  // P[1] = 0
      {3, Rational{7}},  // P[2] = 3
      {1, Rational{8}},  // P[3] = 1
  };

  ASSERT_EQ(result, expected);
}

// ---------------------------------------------------------------------------
// apply_transposed : result = P^T A, row i = A row P[i]
// ---------------------------------------------------------------------------

TEST(PermutationTests, ApplyTransposedPermutesRows) {
  const auto p = reference();
  const Matrix<Rational> expected = {
      {30, 31},  // A row P[0] = 2
      {10, 11},  // A row P[1] = 0
      {40, 41},  // A row P[2] = 3
      {20, 21},  // A row P[3] = 1
  };
  ASSERT_EQ(p.apply_transposed(rows_4x2()), expected);
}

TEST(PermutationTests, ApplyThenApplyTransposedIsIdentity) {
  // P^T (P A) = A
  const auto p = reference();
  const auto A = rows_4x2();
  ASSERT_EQ(p.apply_transposed(p.apply(A)), A);
}

// ---------------------------------------------------------------------------
// post_apply (dense Matrix) : result = A P, column j of A becomes column P[j]
// ---------------------------------------------------------------------------

TEST(PermutationTests, PostApplyDensePermutesColumns) {
  const auto p = reference();
  const Matrix<Rational> A = {
      {0, 1, 2, 3},
      {10, 11, 12, 13},
  };
  const Matrix<Rational> expected = {
      {2, 0, 3, 1},  // columns reordered to P[0..3] = 2,0,3,1
      {12, 10, 13, 11},
  };
  ASSERT_EQ(p.post_apply(A), expected);
}

// ---------------------------------------------------------------------------
// post_apply (scalar) : inverse lookup, returns i such that P[i] == col
// ---------------------------------------------------------------------------

TEST(PermutationTests, PostApplyScalarIsInverse) {
  const auto p = reference();
  ASSERT_EQ(p.post_apply(0), 1u);  // P[1] = 0
  ASSERT_EQ(p.post_apply(1), 3u);  // P[3] = 1
  ASSERT_EQ(p.post_apply(2), 0u);  // P[0] = 2
  ASSERT_EQ(p.post_apply(3), 2u);  // P[2] = 3
}

// ---------------------------------------------------------------------------
// as_dense_matrix : the permutation matrix P with 1 at (P[col], col)
// ---------------------------------------------------------------------------

TEST(PermutationTests, AsDenseMatrix) {
  const auto p = reference();
  const Matrix<Rational> expected = {
      {0, 1, 0, 0},  // 1 at (P[1]=0, 1)
      {0, 0, 0, 1},  // 1 at (P[3]=1, 3)
      {1, 0, 0, 0},  // 1 at (P[0]=2, 0)
      {0, 0, 1, 0},  // 1 at (P[2]=3, 2)
  };
  ASSERT_EQ(p.as_dense_matrix<Rational>(), expected);
}

TEST(PermutationTests, AsDenseMatrixMultiplicationEqualsApply) {
  const auto p = reference();
  const auto A = rows_4x2();
  const Matrix<Rational> dense = p.as_dense_matrix<Rational>();
  ASSERT_EQ(Matrix<Rational>(dense * A), p.apply(A));
}

// ---------------------------------------------------------------------------
// as_sparse_matrix : sparse form of the permutation matrix
// ---------------------------------------------------------------------------

TEST(PermutationTests, AsSparseMatrixMatchesDense) {
  const auto p = reference();
  const Matrix<Rational> from_sparse(p.as_sparse_matrix<Rational>());
  ASSERT_EQ(from_sparse, p.as_dense_matrix<Rational>());
}

// ---------------------------------------------------------------------------
// transposed : the inverse permutation (P^T = P^{-1})
// ---------------------------------------------------------------------------

TEST(PermutationTests, Transposed) {
  const auto p = reference();
  const auto pt = p.transposed();
  // P^{-1}: 0->1, 1->3, 2->0, 3->2
  ASSERT_EQ(pt[0], 1u);
  ASSERT_EQ(pt[1], 3u);
  ASSERT_EQ(pt[2], 0u);
  ASSERT_EQ(pt[3], 2u);
}

TEST(PermutationTests, TransposedAppliedEqualsApplyTransposed) {
  const auto p = reference();
  const auto A = rows_4x2();
  ASSERT_EQ(p.transposed().apply(A), p.apply_transposed(A));
}

TEST(PermutationTests, TransposedTwiceIsOriginal) {
  const auto p = reference();
  const auto ptt = p.transposed().transposed();
  for (size_t i = 0; i < p.size(); ++i) {
    ASSERT_EQ(ptt[i], p[i]);
  }
}

// ---------------------------------------------------------------------------
// is_even : parity of the permutation
// ---------------------------------------------------------------------------

TEST(PermutationTests, IsEvenSingleTransposition) {
  // (1 0) is one transposition -> odd
  ASSERT_FALSE(Permutation::from_vector({1, 0}).is_even());
}

TEST(PermutationTests, IsEvenThreeCycle) {
  // (0 1 2 -> 1 2 0) is a 3-cycle = 2 transpositions -> even
  ASSERT_TRUE(Permutation::from_vector({1, 2, 0}).is_even());
}

TEST(PermutationTests, IsEvenFourCycleIsOdd) {
  // reference() is a single 4-cycle = 3 transpositions -> odd
  ASSERT_FALSE(reference().is_even());
}

// ---------------------------------------------------------------------------
// operator* free functions
// ---------------------------------------------------------------------------

TEST(PermutationTests, OperatorMulLeftEqualsApply) {
  const auto p = reference();
  const auto A = rows_4x2();
  ASSERT_EQ(p * A, p.apply(A));
}

TEST(PermutationTests, OperatorMulRightEqualsPostApply) {
  const auto p = reference();
  const Matrix<Rational> A = {
      {0, 1, 2, 3},
      {10, 11, 12, 13},
  };
  ASSERT_EQ(A * p, p.post_apply(A));
}

// ---------------------------------------------------------------------------
// Randomized tests: random permutation, fixed matrix
// ---------------------------------------------------------------------------

namespace {

// Deterministic fixed matrix of height n and width d, distinct entries.
Matrix<Rational> fixed_matrix(size_t n, size_t d) {
  return Matrix<Rational>::generate(
      n, d, [](size_t i, size_t j) { return Rational(int(i * 100 + j + 1)); });
}

// Builds a Permutation from a shuffled index vector.
Permutation random_permutation(size_t n, std::default_random_engine& rng) {
  std::vector<size_t> indices(n);
  std::iota(indices.begin(), indices.end(), 0);
  std::shuffle(indices.begin(), indices.end(), rng);
  return Permutation::from_vector(std::move(indices));
}

}  // namespace

TEST(PermutationTests, ApplyRandomMatchesDefinition) {
  constexpr size_t n = 9;
  constexpr size_t d = 4;

  const auto A = fixed_matrix(n, d);

  std::default_random_engine rng(12345);

  for (size_t iter = 0; iter < 200; ++iter) {
    SCOPED_TRACE(std::format("iteration: {}", iter));

    const auto p = random_permutation(n, rng);

    // Expected P A: row i of A lands at row P[i] of the result.
    Matrix<Rational> expected(n, d);
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < d; ++j) {
        expected[p[i], j] = A[i, j];
      }
    }

    ASSERT_EQ(p.apply(A), expected);
  }
}

TEST(PermutationTests, ApplyTransposedRandomMatchesDefinition) {
  constexpr size_t n = 9;
  constexpr size_t d = 4;

  const auto A = fixed_matrix(n, d);

  std::default_random_engine rng(67890);

  for (size_t iter = 0; iter < 200; ++iter) {
    SCOPED_TRACE(std::format("iteration: {}", iter));

    const auto p = random_permutation(n, rng);

    // Expected P^T A: row i of the result is row P[i] of A.
    Matrix<Rational> expected(n, d);
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < d; ++j) {
        expected[i, j] = A[p[i], j];
      }
    }

    ASSERT_EQ(p.apply_transposed(A), expected);
  }
}

TEST(PermutationTests, ApplyAndApplyTransposedRoundtripRandom) {
  constexpr size_t n = 9;
  constexpr size_t d = 4;

  const auto A = fixed_matrix(n, d);

  std::default_random_engine rng(13579);

  for (size_t iter = 0; iter < 200; ++iter) {
    SCOPED_TRACE(std::format("iteration: {}", iter));

    const auto p = random_permutation(n, rng);

    // P^T (P A) = A and P (P^T A) = A.
    ASSERT_EQ(p.apply_transposed(p.apply(A)), A);
    ASSERT_EQ(p.apply(p.apply_transposed(A)), A);
  }
}
