#include <gtest/gtest.h>

#include "linalg/Matrix.h"
#include "linalg/expr/TransposedExpr.h"

using namespace linalg;

static_assert(MatrixLike<detail::TransposedExpr<Matrix<double>>, double>);

TEST(TransposedTests, Simple) {
  Matrix<int> matrix = {
      {1, 2},
      {3, 4},
      {5, 6},
  };

  Matrix<int> tr = matrix.transposed();

  Matrix<int> expected = {
      {1, 3, 5},
      {2, 4, 6},
  };

  ASSERT_EQ(tr, expected);
}

TEST(TransposedTests, TransposedConstMatrix) {
  const Matrix<int> matrix = {
      {1, 2},
      {3, 4},
  };

  auto tr = matrix.transposed();
  static_assert(std::same_as<decltype(tr[0, 0]), const int&>);
}
