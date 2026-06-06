#include <gtest/gtest.h>

#include "linalg/CSCMatrix.h"
#include "linalg/Matrix.h"
#include "linalg/expr/SubColsExpr.h"

using namespace linalg;

TEST(SubColsExprTests, Simple) {
  Matrix<int> matrix = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };

  Matrix submatrix = detail::SubColsExpr(matrix, std::vector<size_t>{1});
  Matrix<int> expected = {
      {2},
      {5},
      {8},
  };

  ASSERT_EQ(submatrix, expected);
}
