#include <gtest/gtest.h>

#include "linalg/Assertions.h"
#include "linalg/CSCMatrix.h"
#include "linalg/ConstructSparse.h"
#include "linalg/Matrix.h"

using namespace linalg;

TEST(SumExprTests, DenseSum) {
  Matrix left = {
      {1, 2},
      {3, 4},
  };

  Matrix right = {
      {-5, -10},
      {-20, -30},
  };

  Matrix result = left + right;
  Matrix expected = {
      {-4, -8},
      {-17, -26},
  };

  ASSERT_EQ(result, expected);
}

TEST(SumExprTests, SparseSum) {
  auto left = sparse({
      {0, 1},
      {0, 0},
  });

  auto right = sparse({
      {0, 5},
      {10, 0},
  });

  auto result = left + right;
  std::vector<std::tuple<size_t, size_t, int>> expected = {{0, 1, 6},
                                                           {1, 0, 10}};

  ASSERT_TRIPLES_RANGES_EQ(result.entries(), expected);
}

TEST(SumExprTests, WrongShapes) {
  Matrix left = {
      {1, 2},
  };

  Matrix right = {{3}, {4}};

  ASSERT_ANY_THROW({ left + right; });
}

TEST(SumExprTests, ElementAccess) {
  Matrix left = {
      {1, 2},
      {3, 4},
  };

  Matrix right = {
      {-5, -10},
      {-20, -30},
  };

  auto sum = left + right;

  ASSERT_EQ(sum[0, 0], -4);
  ASSERT_EQ(sum[0, 1], -8);
  ASSERT_EQ(sum[1, 0], -17);
  ASSERT_EQ(sum[1, 1], -26);
}
