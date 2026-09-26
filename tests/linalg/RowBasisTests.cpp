#include <gtest/gtest.h>

#include "../../src/field/BigInteger.h"
#include "linalg/Matrix.h"
#include "linalg/Random.h"
#include "linalg/RowBasis.h"
#include "support/Assertions.h"

using namespace linalg;

TEST(RowBasisTests, UnityMatrix) {
  auto unity = Matrix<Rational>::identity(3);

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

TEST(RowBasisTests, SmallMatrix1) {
  Matrix<Rational> matrix = {
      {-4, 0, -5},
      {-3, -5, 3},
      {-3, -3, 5},
  };

  std::cout << matrix << std::endl;
  auto row_basis = get_row_basis(matrix);

  ASSERT_SETS_EQ(row_basis, (std::vector<size_t>{0, 1, 2}));
}

TEST(RowBasisTests, SmallMatrix2) {
  Matrix<Rational> matrix = {
      {1, 1, 10},
      {2, 2, 20},
  };

  auto row_basis = get_row_basis(matrix);

  ASSERT_EQ(row_basis.size(), 1);
}
