#include <gtest/gtest.h>

#include "linalg/Assertions.h"
#include "linalg/CSCMatrix.h"
#include "linalg/ConstructSparse.h"
#include "linalg/Matrix.h"
#include "linalg/Print.h"
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

TEST(SubColsExprTests, GetColumn) {
  const auto matrix = sparse({
      {1, 2, 3},
      {1, 0, 0},
      {2, 0, 1},
  });

  auto range = detail::SubColsExpr(matrix, std::vector<size_t>{2});

  std::vector<std::pair<size_t, int>> expected = {{0, 3}, {2, 1}};
  ASSERT_DOUBLES_RANGES_EQ(range.col_entries(0), expected);
}

TEST(SubColsExprTests, GetElement) {
  const auto matrix = Matrix<int>{
      {1, 2, 3},
      {1, 0, 0},
      {2, 0, 1},
  };

  auto range = detail::SubColsExpr(matrix, std::vector<size_t>{2});

  ASSERT_EQ((range[0, 0]), 3);
  ASSERT_EQ((range[1, 0]), 0);
  ASSERT_EQ((range[2, 0]), 1);
}


TEST(SubColsExprTests, GetElementManyColumns) {
  const auto matrix = Matrix<int>{
        {1, 2, 3},
        {1, 0, 0},
        {2, 0, 1},
    };

  auto range = detail::SubColsExpr(matrix, std::vector<size_t>{2, 0});

  ASSERT_EQ((range[0, 0]), 3);
  ASSERT_EQ((range[1, 0]), 0);
  ASSERT_EQ((range[2, 0]), 1);

  ASSERT_EQ((range[0, 1]), 1);
  ASSERT_EQ((range[1, 1]), 1);
  ASSERT_EQ((range[2, 1]), 2);
}
