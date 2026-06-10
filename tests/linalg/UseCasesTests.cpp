#include <gtest/gtest.h>

/*
 * Throughout raspisator implementation I use linear algebra a lot. Means must
 * be aligned with the ends not otherwise, so in this test file I compile
 * different tests that highlight the ends not means.
 */

#include <map>

#include "linalg/Linalg.h"

using namespace linalg;

TEST(UseCasesTests, ProblemConstructionVector) {
  Vector<int> cost;

  cost.resize(4);

  for (size_t i = 0; i < 4; ++i) {
    cost[i] = i;
  }

  Vector<int> expected = {0, 1, 2, 3};

  ASSERT_EQ(cost, expected);
}

TEST(UseCasesTests, ProblemConstructionMatrix) {
  CSCMatrix<int> matrix;

  matrix.resize(5, 0);

  for (size_t i = 0; i < 5; ++i) {
    std::map<size_t, int> column;

    column.emplace(0, 321);
    column.emplace(i, 123);

    matrix.add_column(column);
  }

  Matrix<int> expected = {
      {321, 321, 321, 321, 321}, {0, 123, 0, 0, 0}, {0, 0, 123, 0, 0},
      {0, 0, 0, 123, 0},         {0, 0, 0, 0, 123},
  };

  ASSERT_EQ(expected, Matrix(matrix));
}

TEST(UseCasesTests, PermuteCSCRows) {
  CSCMatrix<int> matrix;

  matrix.resize(5, 0);
  matrix.add_column(std::vector<std::pair<size_t, int>>{{0, 1}, {2, 2}});
  matrix.add_column(std::vector<std::pair<size_t, int>>{{1, 3}, {3, 5}});
  matrix.add_column(std::vector<std::pair<size_t, int>>{{4, 6}, {1, 5}});

  // 1 0 0 x
  // 0 3 5 x
  // 2 0 0
  // 0 5 0
  // 0 0 6 x

  std::vector<size_t> rows = {4, 1, 0};

  matrix.map_rows(std::vector<size_t>{2, 1, 4, 4, 0});
  matrix.resize(3, 3);

  ASSERT_EQ(matrix.entries_count(), 4);

  Matrix<int> expected = {
      {0, 0, 6},
      {0, 3, 5},
      {1, 0, 0},
  };

  ASSERT_EQ(expected, Matrix(matrix));
}

TEST(UseCasesTests, GetAdjustedRHS) {
  Vector<int> result(5);

  auto matrix = CSCMatrix<int>::zeros(5);

  matrix.add_column(std::vector<std::pair<size_t, int>>{{0, 1}, {2, 2}});
  matrix.add_column(std::vector<std::pair<size_t, int>>{{1, 3}, {3, 5}});
  matrix.add_column(std::vector<std::pair<size_t, int>>{{4, 6}, {1, 5}});

  result -= matrix.get_column_as_matrix(2) * 5;
  result -= matrix.get_column_as_matrix(1) * -1;

  Vector<int> expected = {0, -22, 0, 5, -30};

  ASSERT_EQ(result, expected);
}

TEST(UseCasesTests, GetBasicCost) {
  Vector<int> cost = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

  std::vector<size_t> basic_vars = {5, 0, 2, 1};

  Vector basic_cost = cost[basic_vars];
  Vector<int> expected = {5, 0, 2, 1};

  ASSERT_EQ(basic_cost, expected);
}

TEST(UseCasesTests, GetReducedCost) {
  // c - A.transposed() * pi

  Vector<int> cost = {-1, -1};

  // 0 2
  // 1 0
  // 3 4
  auto matrix = CSCMatrix<int>::zeros(3);
  matrix.add_column(std::vector<std::pair<size_t, int>>{{1, 1}, {2, 3}});
  matrix.add_column(std::vector<std::pair<size_t, int>>{{0, 2}, {2, 4}});

  Vector<int> pi = {1, 2, 3};

  Vector result = cost - transpose(matrix) * pi;
  Vector<int> expected = {-12, -15};

  ASSERT_EQ(result, expected);
}
