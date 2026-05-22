#include <gtest/gtest.h>

#include "linear/BigInteger.h"
#include "linear/sparse/CSCMatrix.h"

TEST(SparseMatrixTests, ResizeToLarger) {
  CSCMatrix<Rational> matrix = {
      {1, 2, 3},
      {4, 5, 6},
  };

  matrix.resize(5, 5);

  ASSERT_EQ(matrix.rows(), 5);
  ASSERT_EQ(matrix.cols(), 5);

  Matrix<Rational> expected = {
      {1, 2, 3, 0, 0}, {4, 5, 6, 0, 0}, {0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
  };

  ASSERT_EQ(linalg::to_dense(matrix), expected);
}

TEST(SparseMatrixTests, ResizeToSmaller1) {
  CSCMatrix<Rational> matrix = {
      {1, 2, 3},
      {4, 5, 6},
  };

  matrix.resize(1, 1);

  ASSERT_EQ(matrix.rows(), 1);
  ASSERT_EQ(matrix.cols(), 1);

  Matrix<Rational> expected = {{1}};

  ASSERT_EQ(linalg::to_dense(matrix), expected);
}

TEST(SparseMatrixTests, ResizeToSmaller2) {
  CSCMatrix<Rational> matrix = {
      {1, 2, 3},
      {4, 5, 6},
  };

  matrix.resize(1, 3);

  ASSERT_EQ(matrix.rows(), 1);
  ASSERT_EQ(matrix.cols(), 3);

  Matrix<Rational> expected = {{1, 2, 3}};

  ASSERT_EQ(linalg::to_dense(matrix), expected);
}
