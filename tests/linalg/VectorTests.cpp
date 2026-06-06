#include <gtest/gtest.h>

#include "linalg/Vector.h"
#include "linalg/expr/TransposedExpr.h"

using namespace linalg;

static_assert(MatrixLike<Vector<double>, double>);

TEST(VectorTests, Constructor) {
  Vector<int> vector = {1, 2, 3, 4};

  ASSERT_EQ(vector.rows(), 4);
  ASSERT_EQ(vector.cols(), 1);

  ASSERT_EQ(vector.size(), 4);

  for (size_t i = 0; i < 4; ++i) {
    ASSERT_EQ(vector[i], i + 1);
    ASSERT_EQ((vector[i, 0]), i + 1);
  }
}

TEST(VectorTests, AssignTransposedVectorToMatrix) {
  Vector<int> vector = {1, 2, 3, 4};

  Matrix<int> matrix = vector.transposed();
  Matrix<int> expected = {{1, 2, 3, 4}};

  ASSERT_EQ(matrix, expected);
}
