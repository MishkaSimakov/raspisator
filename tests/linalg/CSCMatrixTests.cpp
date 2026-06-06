#include <gtest/gtest.h>

#include "linalg/CSCMatrix.h"
#include "linalg/Matrix.h"
#include "linalg/Print.h"

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

TEST(CSCMatrixTests, ConstructByColumnsDuplicatedRow) {
  auto matrix = CSCMatrix<int>::zeros(3, 0);

  matrix.add_column(
      std::vector<std::pair<size_t, int>>{{0, 1}, {2, 2}, {0, 5}, {0, 6}});

  Matrix<int> dense = matrix;
  Matrix<int> expected = {
      {12},
      {0},
      {2},
  };

  std::cout << dense << std::endl;

  ASSERT_EQ(dense, expected);
}

TEST(CSCMatrixTests, ConstructByColumnsRowOutside) {
  auto matrix = CSCMatrix<int>::zeros(3, 0);

  ASSERT_ANY_THROW(
      { matrix.add_column(std::vector<std::pair<size_t, int>>{{5, 5}}); });
}

TEST(CSCMatrixTests, Subscript) {
  auto matrix = CSCMatrix<int>::zeros(5, 0);

  matrix.add_column(std::vector<std::pair<size_t, int>>{{0, 1}, {2, 3}});
  matrix.add_column(std::vector<std::pair<size_t, int>>{{1, 2}, {3, 4}});

  Matrix<int> expected = {
      {1, 0}, {0, 2}, {3, 0}, {0, 4}, {0, 0},
  };

  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < 2; ++j) {
      ASSERT_EQ((matrix[i, j]), (expected[i, j]));
    }
  }
}

TEST(CSCMatrixTests, Transposed) {
  auto matrix = CSCMatrix<int>::zeros(5, 0);

  matrix.add_column(std::vector<std::pair<size_t, int>>{{0, 1}, {2, 3}});
  matrix.add_column(std::vector<std::pair<size_t, int>>{{1, 2}, {3, 4}});

  Matrix dense = matrix.transposed();

  Matrix<int> expected = {
      {1, 0, 3, 0, 0},
      {0, 2, 0, 4, 0},
  };

  ASSERT_EQ(dense, expected);
}
