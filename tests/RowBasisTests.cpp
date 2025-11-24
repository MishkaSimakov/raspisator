#include <gtest/gtest.h>

#include "../src/linear/matrix/RowBasis.h"
#include "Assertions.h"
#include "linear/BigInteger.h"
#include "linear/matrix/Matrix.h"

TEST(RowBasisTests, UnityMatrix) {
  auto unity = Matrix<Rational>::unity(3);

  auto row_basis = linalg::get_row_basis(unity);

  ASSERT_SETS_EQ(row_basis, (std::vector<size_t>{0, 1, 2}));
}

TEST(RowBasisTests, WideMatrix) {
  auto matrix = Matrix<Rational>(3, 4, 0);
  for (size_t i = 0; i < 3; ++i) {
    matrix[i, i] = 1;
  }

  auto row_basis = linalg::get_row_basis(matrix);

  ASSERT_SETS_EQ(row_basis, (std::vector<size_t>{0, 1, 2}));
}

TEST(RowBasisTests, LinearlyDependentRows) {
  Matrix<Rational> matrix = {{1, 1, 0, 0}, {0, 0, 1, 1}, {1, 1, 1, 1}};

  auto row_basis = linalg::get_row_basis(matrix);

  ASSERT_SETS_EQ(row_basis, (std::vector<size_t>{0, 1}));
}

TEST(RowBasisTests, ZeroRow) {
  Matrix<Rational> matrix = {{1, 0}, {0, 1}, {0, 0}};

  auto row_basis = linalg::get_row_basis(matrix);

  ASSERT_SETS_EQ(row_basis, (std::vector<size_t>{0, 1}));
}

TEST(RowBasisTests, LongRows) {
  Matrix<Rational> matrix = {
      {1, 1, 0, 0}, {0, 0, 1, 1}, {0, 1, 0, 1}, {1, 0, 1, 0}};

  auto row_basis = linalg::get_row_basis(matrix);

  ASSERT_SETS_EQ(row_basis, (std::vector<size_t>{0, 1, 2}));
}

TEST(RowBasisTests, RowPermutations) {
  Matrix<Rational> matrix = {
      {1, 0},
      {0, 0},
      {0, 1},
  };

  auto row_basis = linalg::get_row_basis(matrix);

  ASSERT_SETS_EQ(row_basis, (std::vector<size_t>{0, 2}));
}

TEST(RowBasisTests, CompleteRowBasis1) {
  Matrix<Rational> matrix = {
      {1, 0},
      {0, 0},
      {0, 1},
  };

  auto row_basis = linalg::complete_row_basis(matrix, std::vector<size_t>{0});

  ASSERT_SETS_EQ(row_basis, (std::vector<size_t>{0, 2}));
}

TEST(RowBasisTests, CompleteRowBasis2) {
  Matrix<Rational> matrix = {
    {1, 0},
    {0, 0},
    {0, 1},
    {0, 0}
};

  auto row_basis = linalg::complete_row_basis(matrix, std::vector<size_t>{0});

  ASSERT_SETS_EQ(row_basis, (std::vector<size_t>{0, 2}));
}

TEST(RowBasisTests, CompleteRowBasisUnityMatrix) {
  auto matrix = Matrix<Rational>::unity(50);

  std::vector<size_t> expected(50);
  std::iota(expected.begin(), expected.end(), 0);

  {
    std::vector<size_t> partial_basis(25);
    std::iota(partial_basis.begin(), partial_basis.end(), 0);

    auto row_basis = linalg::complete_row_basis(matrix, partial_basis);

    ASSERT_SETS_EQ(row_basis, expected);
  }

  {
    auto row_basis = linalg::complete_row_basis(matrix, expected);
    ASSERT_SETS_EQ(row_basis, expected);
  }
}
