#include <gtest/gtest.h>

#include "linear/BigInteger.h"
#include "linear/Matrix.h"
#include "linear/RowBasis.h"

void assert_sets_equal(auto&& left_range, auto&& right_range) {
  auto left_set =
      std::set(std::from_range, std::forward<decltype(left_range)>(left_range));

  auto right_set = std::set(std::from_range,
                            std::forward<decltype(right_range)>(right_range));

  ASSERT_EQ(left_set, right_set);
}

#define ASSERT_SETS_EQ(left, right) assert_sets_equal(left, right)

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
