#include <gtest/gtest.h>

#include "linalg/Matrix.h"
#include "linalg/expr/TransposedExpr.h"

using namespace linalg;

static_assert(MatrixLike<TransposedExpr<Matrix<double>>, double>);

TEST(TransposedTests, Simple) {
  Matrix<int> matrix = {
      {1, 2},
      {3, 4},
      {5, 6},
  };

  Matrix<int> tr = transposed(matrix);

  Matrix<int> expected = {
      {1, 3, 5},
      {2, 4, 6},
  };

  ASSERT_EQ(tr, expected);
}

TEST(TransposedTests, Assignment) {
  Matrix<int> matrix = {
      {1, 2},
      {3, 4},
  };

  auto tr = transposed(matrix);
  tr[1, 0] = 123;

  Matrix<int> expected = {
      {1, 123},
      {3, 4},
  };

  ASSERT_EQ(matrix, expected);
}

TEST(TransposedTests, TransposedConstMatrix) {
  const Matrix<int> matrix = {
      {1, 2},
      {3, 4},
  };

  auto tr = transposed(matrix);
  static_assert(std::same_as<decltype(tr[0, 0]), const int&>);
}
