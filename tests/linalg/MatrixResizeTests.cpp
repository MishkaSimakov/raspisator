#include <gtest/gtest.h>

#include <random>

#include "linalg/Matrix.h"

using namespace linalg;

TEST(MatrixResizeTests, Extend) {
  Matrix<int> matrix = {{1}};

  matrix.resize(2, 2);

  Matrix<int> expected = {
      {1, 0},
      {0, 0},
  };

  ASSERT_EQ(matrix, expected);
}

TEST(MatrixResizeTests, Shrink) {
  Matrix<int> matrix = {
      {1, 2},
      {3, 4},
  };

  matrix.resize(1, 1);

  Matrix<int> expected = {{1}};

  ASSERT_EQ(matrix, expected);
}

TEST(MatrixResizeTests, Random) {
  std::default_random_engine random;

  std::uniform_int_distribution<size_t> random_size(0, 25);
  std::uniform_int_distribution<int> random_value(-100, 100);

  for (size_t i = 0; i < 10'000; ++i) {
    const size_t old_rows = random_size(random);
    const size_t old_cols = random_size(random);
    const size_t new_rows = random_size(random);
    const size_t new_cols = random_size(random);

    auto matrix = Matrix<int>::generate(
        old_rows, old_cols,
        [&](size_t, size_t) { return random_value(random); });

    auto copy = matrix;

    matrix.resize(new_rows, new_cols);

    ASSERT_EQ(matrix.rows(), new_rows);
    ASSERT_EQ(matrix.cols(), new_cols);

    for (size_t row = 0; row < new_rows; ++row) {
      for (size_t col = 0; col < new_cols; ++col) {
        if (row < old_rows && col < old_cols) {
          ASSERT_EQ((matrix[row, col]), (copy[row, col]));
        } else {
          ASSERT_EQ((matrix[row, col]), 0);
        }
      }
    }
  }
}
