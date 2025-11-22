#include <gtest/gtest.h>

#include "linear/BigInteger.h"
#include "linear/matrix/Matrix.h"

// some type checks
Matrix<double> get_matrix() { std::unreachable(); }
const Matrix<double>& get_const_matrix() { std::unreachable(); }

static_assert(std::same_as<decltype(get_matrix()[0, 0]), double&>);
static_assert(std::same_as<decltype(get_const_matrix()[0, 0]), const double&>);

static_assert(std::same_as<decltype(get_matrix()[std::pair{0, 1}, 1]),
                           MatrixSlice<double>>);
static_assert(std::same_as<decltype(get_const_matrix()[std::pair{0, 1}, 0]),
                           MatrixSlice<const double>>);

static_assert(std::same_as<decltype(get_matrix()[0, std::pair{0, 1}]),
                           MatrixSlice<double>>);
static_assert(std::same_as<decltype(get_const_matrix()[0, std::pair{0, 1}]),
                           MatrixSlice<const double>>);

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
  Matrix<Rational> matrix = {{5}};
  matrix.inverse();

  Matrix<Rational> expected = {{Rational(1) / 5}};

  ASSERT_EQ(matrix, expected);
}

TEST(MatrixTests, InverseMatrix) {
  Matrix<Rational> matrix = {{4, 3}, {3, 2}};
  matrix.inverse();

  Matrix<Rational> expected = {{-2, 3}, {3, -4}};

  ASSERT_EQ(matrix, expected);
}

TEST(MatrixTests, InvertMultiply) {
  {
    Matrix<Rational> matrix = {{4, 3}, {3, 2}};

    auto inverse = matrix;
    inverse.inverse();

    ASSERT_EQ(matrix * inverse, Matrix<Rational>::unity(2));
    ASSERT_EQ(inverse * matrix, Matrix<Rational>::unity(2));
  }

  {
    Matrix<Rational> matrix = {{0, 3, 0}, {0, 0, 2}, {1, 0, 0}};

    auto inverse = matrix;
    inverse.inverse();

    ASSERT_EQ(matrix * inverse, Matrix<Rational>::unity(3));
    ASSERT_EQ(inverse * matrix, Matrix<Rational>::unity(3));
  }
}

TEST(MatrixTests, AddSubtract) {
  Matrix<int> first = {{4, 3}, {3, 2}};
  Matrix<int> second = {{5, 6}, {7, 8}};

  ASSERT_EQ(first + second - second, first);
  ASSERT_EQ(second + first - first, second);
  ASSERT_EQ(first - second + second, first);
}

TEST(MatrixTests, SimpleGaussElimination) {
  Matrix<Rational> matrix = {{1, 0}, {1, 1}};

  matrix.gaussian_elimination(0, 0);

  ASSERT_EQ(matrix, (Matrix<Rational>{{1, 0}, {0, 1}}));
}

TEST(MatrixTests, GaussEliminationWithZeros) {
  Matrix<Rational> matrix = {{1, 0, 1}, {1, 1, 1}, {0, 1, 1}};

  matrix.gaussian_elimination(0, 0);

  ASSERT_EQ(matrix, (Matrix<Rational>{{1, 0, 1}, {0, 1, 0}, {0, 1, 1}}));
}

TEST(MatrixTests, vstack) {
  Matrix top = {{1, 2}};
  Matrix bottom = {{3, 4}, {5, 6}};

  {
    Matrix expected = {{1, 2}, {3, 4}, {5, 6}};
    ASSERT_EQ(linalg::vstack(top, bottom), expected);
  }

  {
    Matrix expected = {{1, 2}, {3, 4}};
    ASSERT_EQ(linalg::vstack(top, bottom[0, {0, 2}]), expected);
  }

  {
    Matrix expected = {{1, 2}, {5, 6}};
    ASSERT_EQ(linalg::vstack(top, bottom[1, {0, 2}]), expected);
  }
}

TEST(MatrixTests, hstack) {
  Matrix left = {{1, 2}, {3, 4}};
  Matrix right = {{5, 6}, {7, 8}};

  {
    Matrix expected = {{1, 2, 5, 6}, {3, 4, 7, 8}};
    ASSERT_EQ(linalg::hstack(left, right), expected);
  }

  {
    Matrix expected = {{1, 5}, {3, 7}};
    ASSERT_EQ(linalg::hstack(left[{0, 2}, 0], right[{0, 2}, 0]), expected);
  }
}

TEST(MatrixTests, SliceAddition) {
  Matrix left = {{1, 2}, {3, 4}};
  Matrix right = {{5, 6}, {7, 8}};

  left[0, {0, 2}] += right[0, {0, 2}];
  ASSERT_EQ(left, (Matrix{{6, 8}, {3, 4}}));

  left[0, {0, 2}] -= left[0, {0, 2}];
  ASSERT_EQ(left, (Matrix{{0, 0}, {3, 4}}));

  left[0, {0, 2}] += left[1, {0, 2}];
  ASSERT_EQ(left, (Matrix{{3, 4}, {3, 4}}));
}

TEST(MatrixTests, SliceMultiplication) {
  Matrix matrix = {{1, 2}, {3, 4}};

  matrix[{0, 2}, 1] *= 42;

  ASSERT_EQ(matrix, (Matrix{{1, 84}, {3, 168}}));
}

TEST(MatrixTests, MatrixFromSlice) {
  Matrix A = {{1, 2}, {3, 4}};

  {
    Matrix sliced = A[{0, 2}, 0];
    ASSERT_EQ(sliced, (Matrix{{1}, {3}}));
  }

  {
    Matrix sliced = A[1, {0, 2}];
    ASSERT_EQ(sliced, (Matrix{{3, 4}}));
  }

  {
    const Matrix<int>& const_ref = A;
    Matrix sliced = const_ref[1, {0, 2}];
    ASSERT_EQ(sliced, (Matrix{{3, 4}}));
  }
}

TEST(MatrixTests, SliceFromSlice) {
  Matrix A = {{1, 2}, {3, 4}};

  auto slice = A[{0, 2}, 1];
  ASSERT_EQ((slice[0, 0]), 2);
  ASSERT_EQ((slice[1, 0]), 4);
}

TEST(MatrixTests, SliceFromSliceAssignment) {
  Matrix A = {{1, 2}, {3, 4}};

  auto slice = A[{0, 2}, 1];
  slice[0, 0] = 123;

  ASSERT_EQ((A[0, 1]), 123);
}

TEST(MatrixTests, SliceAssignment) {
  Matrix A = {{0}};
  Matrix B = {{0, 123}};

  auto B_slice = B[{0, 1}, 1];

  A[0, {0, 1}] = B_slice[0, {0, 1}];

  ASSERT_EQ(A, (Matrix{{123}}));
}
