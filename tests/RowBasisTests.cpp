#include <gtest/gtest.h>

#include "Assertions.h"
#include "linear/BigInteger.h"
#include "linear/Matrix.h"
#include "linear/RowBasis.h"

TEST(RowBasisTests, UnityMatrix) {
  auto unity = Matrix<Rational>::unity(3);

  auto row_basis = get_row_basis(unity);

  ASSERT_SETS_EQ(row_basis, (std::vector<size_t>{0, 1, 2}));
}

TEST(RowBasisTests, LinearlyDependentRows) {
  Matrix<Rational> matrix = {{1, 1, 0, 0}, {0, 0, 1, 1}, {1, 1, 1, 1}};

  auto row_basis = get_row_basis(matrix);

  ASSERT_SETS_EQ(row_basis, (std::vector<size_t>{0, 1}));
}

TEST(RowBasisTests, ZeroRow) {
  Matrix<Rational> matrix = {{1, 0}, {0, 1}, {0, 0}};

  auto row_basis = get_row_basis(matrix);

  ASSERT_SETS_EQ(row_basis, (std::vector<size_t>{0, 1}));
}

TEST(RowBasisTests, LongRows) {
  Matrix<Rational> matrix = {
      {1, 1, 0, 0}, {0, 0, 1, 1}, {0, 1, 0, 1}, {1, 0, 1, 0}};

  auto row_basis = get_row_basis(matrix);

  ASSERT_SETS_EQ(row_basis, (std::vector<size_t>{0, 1, 2}));
}

TEST(RowBasisTests, RowPermutations) {
  Matrix<Rational> matrix = {
      {1, 0},
      {0, 0},
      {0, 1},
  };

  auto row_basis = get_row_basis(matrix);

  ASSERT_SETS_EQ(row_basis, (std::vector<size_t>{0, 2}));
}
