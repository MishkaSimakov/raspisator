#include <gtest/gtest.h>

#include "TestMatrices.h"
#include "linear/BigInteger.h"
#include "linear/matrix/Elimination.h"
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

  ASSERT_NO_FATAL_FAILURE(check_L(L));
  ASSERT_NO_FATAL_FAILURE(check_U(U));

  ASSERT_EQ(L * U, matrix);
}

TEST_P(CommonLUTests, SparseDecomposeThenCompose) {
  auto [_, matrix] = GetParam();

  auto [n, d] = matrix.shape();
  auto sparse = CSCMatrix<Rational>(matrix);

  std::vector<size_t> columns(matrix.get_width());
  std::iota(columns.begin(), columns.end(), 0);

  auto [P, Q, ls, us] =
      linalg::FullPivotingLU<Rational>(n).get(sparse, columns);

  auto L = Matrix<Rational>::unity(n);
  for (auto entry : ls) {
    L = ls.apply(std::move(L), entry);
  }

  auto U = Matrix<Rational>::unity(n);
  for (auto entry : us | std::views::reverse) {
    U = us.apply(std::move(U), entry);
  }

  ASSERT_NO_FATAL_FAILURE(check_L(L));
  ASSERT_NO_FATAL_FAILURE(check_U(U));

  auto expected = linalg::to_dense(sparse);

  ASSERT_EQ(Q * U * L * P * expected, Matrix<Rational>::unity(n));
}

INSTANTIATE_TEST_SUITE_P(
    CommonLUTests, CommonLUTests, ::testing::ValuesIn(test_matrices()),
    [](const testing::TestParamInfo<CommonLUTests::ParamType>& info) {
      return info.param.first;
    });
