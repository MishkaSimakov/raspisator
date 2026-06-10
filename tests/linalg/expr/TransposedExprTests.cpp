#include <gtest/gtest.h>

#include "../../ConstructSparse.h"
#include "linalg/Matrix.h"
#include "linalg/MatrixAssertions.h"
#include "linalg/MatrixSpy.h"
#include "linalg/Transpose.h"
#include "linalg/expr/TransposedExpr.h"

using namespace linalg;

static_assert(MatrixRange<detail::TransposedExpr<Matrix<double>>>);

TEST(TransposedExprTests, Simple) {
  Matrix<int> matrix = {
      {1, 2},
      {3, 4},
      {5, 6},
  };

  Matrix<int> tr = transpose(matrix);

  Matrix<int> expected = {
      {1, 3, 5},
      {2, 4, 6},
  };

  ASSERT_EQ(tr, expected);
}

TEST(TransposedExprTests, ElementWiseAccess) {
  Matrix<int> matrix = {
      {1, 2},
      {3, 4},
      {5, 6},
  };

  auto tr = transpose(matrix);

  ASSERT_EQ((tr[0, 0]), 1);
  ASSERT_EQ((tr[0, 1]), 3);
  ASSERT_EQ((tr[0, 2]), 5);
  ASSERT_EQ((tr[1, 0]), 2);
  ASSERT_EQ((tr[1, 1]), 4);
  ASSERT_EQ((tr[1, 2]), 6);
}

TEST(TransposedExprTests, ColWiseMatrix) {
  auto matrix = sparse({
      {0, 1, 2},
      {0, 0, 3},
  });

  auto tr = transpose(matrix);

  static_assert(RowWiseMatrixRange<decltype(tr)>);

  std::vector<std::pair<size_t, int>> expected_row1 = {{0, 1}};
  ASSERT_ROW_ENTRIES_EQ(tr, 1, expected_row1);

  std::vector<std::pair<size_t, int>> expected_row2 = {{0, 2}, {1, 3}};
  ASSERT_ROW_ENTRIES_EQ(tr, 2, expected_row2);
}

TEST(TransposedExprTests, RowWiseMatrixRange) {
  auto matrix = sparse({
      {0, 1, 2},
      {0, 0, 3},
  });

  auto tr = detail::TransposedExpr<
      detail::TransposedExpr<detail::RefView<CSCMatrix<int>>>>(
      transpose(matrix));

  static_assert(ColWiseMatrixRange<decltype(tr)>);

  std::vector<std::pair<size_t, int>> expected_col1 = {{0, 1}};
  ASSERT_COL_ENTRIES_EQ(tr, 1, expected_col1);

  std::vector<std::pair<size_t, int>> expected_col2 = {{0, 2}, {1, 3}};
  ASSERT_COL_ENTRIES_EQ(tr, 2, expected_col2);
}

TEST(TransposedExprTests, TransposeFunctionDoesntCopy) {
  MatrixSpy spy;

  MatrixSpy::reset_counters();
  auto tr_view = transpose(detail::TransposedExpr(std::move(spy)));

  ASSERT_EQ(MatrixSpy::copy_constructor_calls, 0);
  ASSERT_EQ(MatrixSpy::copy_assignment_calls, 0);
}

TEST(TransposedExprTests, DoubleTransposeFunctionWithMatrix) {
  Matrix matrix = {
      {1, 2, 3},
      {4, 5, 6},
  };

  auto result = transpose(transpose(matrix));

  ASSERT_EQ(result, matrix);
}
