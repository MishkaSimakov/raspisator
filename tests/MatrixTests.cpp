#include <gtest/gtest.h>

#include <iostream>

#include "linear/Matrix.h"

TEST(MatrixTests, SimpleInitialization) {
  Matrix<int> matrix(2, 3);

  auto [n, m] = matrix.shape();
  ASSERT_EQ(n, 2);
  ASSERT_EQ(m, 3);

  ASSERT_EQ((matrix[0, 0]), 0);
  ASSERT_EQ((matrix[1, 0]), 0);
}

TEST(MatrixTests, ElementsAssignment) {
  Matrix<int> matrix(2, 3);

  matrix[1, 2] = 123;

  ASSERT_EQ((matrix[1, 2]), 123);
}

TEST(MatrixTests, InitializerList) {
  Matrix<int> matrix1 = {{1, 2}};

  ASSERT_EQ((matrix1[0, 0]), 1);
  ASSERT_EQ((matrix1[0, 1]), 2);

  Matrix<int> matrix2 = {{1}, {2}};

  ASSERT_EQ((matrix2[0, 0]), 1);
  ASSERT_EQ((matrix2[1, 0]), 2);
}

TEST(MatrixTests, Equality) {
  Matrix<int> matrix1 = {{1, 2}};
  Matrix<int> matrix2 = {{1, 2}};
  Matrix<int> matrix3 = {{1, 2}, {3, 4}};
  Matrix<int> matrix4 = {{2, 2}};

  ASSERT_EQ(matrix1, matrix2);
  ASSERT_NE(matrix1, matrix3);
  ASSERT_NE(matrix1, matrix4);
}

TEST(MatrixTests, InverseScalar) {
  Matrix<double> matrix = {{5}};
  matrix.inverse();

  Matrix<double> expected = {{0.2}};

  ASSERT_EQ(matrix, expected);
}

TEST(MatrixTests, InverseMatrix) {
  Matrix<double> matrix = {{4, 3}, {3, 2}};
  matrix.inverse();

  Matrix<double> expected = {{-2, 3}, {3, -4}};

  ASSERT_EQ(matrix, expected);
}

TEST(MatrixTests, InvertMultiply) {
  Matrix<double> matrix = {{4, 3}, {3, 2}};

  auto inverse = matrix;
  inverse.inverse();

  ASSERT_EQ(matrix * inverse, Matrix<double>::unity(2));
  ASSERT_EQ(inverse * matrix, Matrix<double>::unity(2));
}

TEST(MatrixTests, AddSubtract) {
  Matrix<int> first = {{4, 3}, {3, 2}};
  Matrix<int> second = {{5, 6}, {7, 8}};

  ASSERT_EQ(first + second - second, first);
  ASSERT_EQ(second + first - first, second);
  ASSERT_EQ(first - second + second, first);
}
