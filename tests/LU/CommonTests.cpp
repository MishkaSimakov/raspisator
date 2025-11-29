#include <gtest/gtest.h>

#include "TestMatrices.h"
#include "linear/BigInteger.h"
#include "linear/matrix/LU.h"
#include "linear/matrix/Matrix.h"
#include "linear/sparse/CSCMatrix.h"
#include "linear/sparse/LU.h"

template <typename Field>
void check_U(const Matrix<Field>& matrix) {
  auto [n, d] = matrix.shape();

  ASSERT_EQ(n, d);

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < i; ++j) {
      ASSERT_EQ((matrix[i, j]), 0);
    }
  }
}

template <typename Field>
void check_L(const Matrix<Field>& matrix) {
  auto [n, d] = matrix.shape();

  ASSERT_EQ(n, d);

  for (size_t i = 0; i < n; ++i) {
    ASSERT_EQ((matrix[i, i]), 1);

    for (size_t j = i + 1; j < n; ++j) {
      ASSERT_EQ((matrix[i, j]), 0);
    }
  }
}

class CommonLUTests : public ::testing::TestWithParam<
                          std::pair<std::string, Matrix<Rational>>> {
 protected:
  Matrix<Rational> matrix_;
};

TEST_P(CommonLUTests, DenseDecomposeThenCompose) {
  auto [_, matrix] = GetParam();

  auto [L, U] = linalg::get_lu(matrix);

  check_L(L);
  check_U(U);

  ASSERT_EQ(L * U, matrix);
}

TEST_P(CommonLUTests, SparseDecomposeThenCompose) {
  auto [_, matrix] = GetParam();
  auto sparse = CSCMatrix<Rational>(matrix);

  auto [L, U, P] = linalg::sparse_lu(sparse);
  auto dense_L = linalg::to_dense(L);
  auto dense_U = linalg::to_dense(U);

  // add ones on the diagonal
  for (size_t i = 0; i < matrix.get_height(); ++i) {
    dense_L[i, i] = 1;
  }

  check_L(dense_L);
  check_U(dense_U);

  auto expected = linalg::to_dense(linalg::apply_permutation(sparse, P));

  ASSERT_EQ(dense_L * dense_U, expected);
}

INSTANTIATE_TEST_SUITE_P(
    CommonLUTests, CommonLUTests, ::testing::ValuesIn(test_matrices()),
    [](const testing::TestParamInfo<CommonLUTests::ParamType>& info) {
      return info.param.first;
    });
