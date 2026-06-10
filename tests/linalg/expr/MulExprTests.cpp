#include <gtest/gtest.h>

#include "ConstructSparse.h"
#include "linalg/Linalg.h"
#include "linalg/MatrixSpy.h"

using namespace linalg;

TEST(MulExprTests, DenseMatrixDenseMatrix) {
  Matrix left = {
      {1, 2, 3},
      {4, 5, 6},
  };

  Matrix right = {
      {10, 20, 30, 40},
      {50, 60, 70, 80},
      {100, 200, 300, 400},
  };

  Matrix expected = {
      {410, 740, 1070, 1400},
      {890, 1580, 2270, 2960},
  };

  Matrix result = left * right;

  ASSERT_EQ(result, expected);
}

TEST(MulExprTests, DenseMatrixDenseVector) {
  Matrix left = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };

  Vector right = {-5, 6, 100};

  Vector expected = {307, 610, 913};

  Matrix result = left * right;

  ASSERT_EQ(result, expected);
}

TEST(MulExprTests, DenseVectorDenseMatrix) {
  Vector left = {-5, 6, 100};

  Matrix right = {
      {1, 4, 7},
      {2, 5, 8},
      {3, 6, 9},
  };

  Matrix expected = {{307, 610, 913}};

  Matrix result = transpose(left) * right;

  ASSERT_EQ(result, expected);
}

TEST(MulExprTests, ColSparseDenseMatrix) {
  auto left = sparse({
      {-1, 0, 0},
      {0, 2, 1},
      {3, 0, 1},
      {4, 1, 2},
  });

  Matrix right = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };

  Matrix expected = {
      {-1, -2, -3},
      {15, 18, 21},
      {10, 14, 18},
      {22, 29, 36},
  };

  Matrix result = left * right;

  ASSERT_EQ(result, expected);
}

TEST(MulExprTests, DenseMatrixRowSparse) {
  Matrix left = {
      {1, 4, 7},
      {2, 5, 8},
      {3, 6, 9},
  };

  auto right = sparse({
      {-1, 0, 0},
      {0, 2, 1},
      {3, 0, 1},
      {4, 1, 2},
  });

  Matrix expected = {
      {-1, 15, 10, 22},
      {-2, 18, 14, 29},
      {-3, 21, 18, 36},
  };

  Matrix result = left * transpose(right);

  ASSERT_EQ(result, expected);
}

TEST(MulExprTests, DenseMatrixColSparse) {
  Matrix left = {
      {1, 4, 7},
      {2, 5, 8},
      {3, 6, 9},
  };

  auto right = sparse({
      {-1, 0, 3, 4},
      {0, 2, 0, 1},
      {0, 1, 1, 2},
  });

  Matrix expected = {
      {-1, 15, 10, 22},
      {-2, 18, 14, 29},
      {-3, 21, 18, 36},
  };

  Matrix result = left * right;

  ASSERT_EQ(result, expected);
}

TEST(MulExprTests, DoNotCopyMatrix) {
  MatrixSpy left;
  MatrixSpy right;

  MatrixSpy::reset_counters();

  auto view = left * right;

  ASSERT_EQ(MatrixSpy::copy_constructor_calls, 0);
  ASSERT_EQ(MatrixSpy::move_constructor_calls, 0);
  ASSERT_EQ(MatrixSpy::copy_assignment_calls, 0);
  ASSERT_EQ(MatrixSpy::move_assignment_calls, 0);
}

TEST(MulExprTests, WrapsTemporariesInOwningView) {
  MatrixSpy left;
  MatrixSpy right;

  auto view1 = left * MatrixSpy();
  static_assert(std::same_as<typename decltype(view1)::RightType,
                             detail::OwningView<MatrixSpy>>);

  auto view2 = MatrixSpy() * right;
  static_assert(std::same_as<typename decltype(view2)::LeftType,
                             detail::OwningView<MatrixSpy>>);
}

TEST(MulExprTests, WrongShapes) {
  Matrix left = {
      {1, 2},
  };

  Matrix right = {
      {1, 2},
  };

  ASSERT_ANY_THROW({ left* right; });
}
