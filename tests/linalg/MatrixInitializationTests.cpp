#include <gtest/gtest.h>

#include "linalg/Matrix.h"

using namespace linalg;

TEST(MatrixInitializationTests, Zeros) {
  auto matrix = Matrix<int>::zeros(3, 3);

  ASSERT_EQ(matrix.rows(), 3);
  ASSERT_EQ(matrix.cols(), 3);

  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      ASSERT_EQ((matrix[i, j]), 0);
    }
  }
}

TEST(MatrixInitializationTests, Uninitialized) {
  auto matrix = Matrix<int>::uninitialized(3, 3);

  ASSERT_EQ(matrix.rows(), 3);
  ASSERT_EQ(matrix.cols(), 3);
}

TEST(MatrixInitializationTests, InitializerList) {
  Matrix<int> matrix = {
      {1, 2},
      {3, 4},
  };

  ASSERT_EQ(matrix.rows(), 2);
  ASSERT_EQ(matrix.cols(), 2);

  ASSERT_EQ((matrix[0, 0]), 1);
  ASSERT_EQ((matrix[0, 1]), 2);
  ASSERT_EQ((matrix[1, 0]), 3);
  ASSERT_EQ((matrix[1, 1]), 4);
}

TEST(MatrixInitializationTests, WrongInitializerListShape) {
  ASSERT_ANY_THROW(({ Matrix<int> matrix = {{1, 2, 3}, {1, 2}}; }));
}

TEST(MatrixInitializationTests, Generate) {
  auto generator = [](size_t i, size_t j) {
    return 10 * i + j;
  };

  auto matrix = Matrix<int>::generate(5, 5, generator);

  ASSERT_EQ(matrix.rows(), 5);
  ASSERT_EQ(matrix.cols(), 5);

  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < 5; ++j) {
      ASSERT_EQ((matrix[i, j]), generator(i, j));
    }
  }
}
