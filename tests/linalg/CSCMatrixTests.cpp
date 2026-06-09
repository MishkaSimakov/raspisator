#include <gtest/gtest.h>

#include "linalg/CSCMatrix.h"
#include "linalg/Matrix.h"
#include "linalg/Print.h"
#include "linalg/Transpose.h"

using namespace linalg;

TEST(CSCMatrixTests, ConstructByColumns) {
  auto matrix = CSCMatrix<int>::zeros(5, 0);

  matrix.add_column(std::vector<std::pair<size_t, int>>{{0, 1}, {2, 3}});
  matrix.add_column(std::vector<std::pair<size_t, int>>{{1, 2}, {3, 4}});

  Matrix<int> dense = matrix;
  Matrix<int> expected = {
      {1, 0}, {0, 2}, {3, 0}, {0, 4}, {0, 0},
  };

  ASSERT_EQ(dense, expected);
}

TEST(CSCMatrixTests, ConstructByColumnsRowOutside) {
  auto matrix = CSCMatrix<int>::zeros(3, 0);

  ASSERT_ANY_THROW(
      { matrix.add_column(std::vector<std::pair<size_t, int>>{{5, 5}}); });
}

TEST(CSCMatrixTests, Transposed) {
  auto matrix = CSCMatrix<int>::zeros(5, 0);

  matrix.add_column(std::vector<std::pair<size_t, int>>{{0, 1}, {2, 3}});
  matrix.add_column(std::vector<std::pair<size_t, int>>{{1, 2}, {3, 4}});

  Matrix dense = transpose(matrix);

  Matrix<int> expected = {
      {1, 0, 3, 0, 0},
      {0, 2, 0, 4, 0},
  };

  ASSERT_EQ(dense, expected);
}
