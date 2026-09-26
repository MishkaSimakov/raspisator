#include <gtest/gtest.h>

#include "ConstructSparse.h"
#include "linalg/Linalg.h"
#include "support/GMPRational.h"

using namespace linalg;

TEST(SparseMatrixTests, ResizeToLarger) {
  auto matrix = sparse<GMPRational>({
      {1, 2, 3},
      {4, 5, 6},
  });

  matrix.resize(5, 5);

  ASSERT_EQ(matrix.rows(), 5);
  ASSERT_EQ(matrix.cols(), 5);

  Matrix<GMPRational> expected = {
      {1, 2, 3, 0, 0}, {4, 5, 6, 0, 0}, {0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
  };

  ASSERT_EQ(Matrix(matrix), expected);
}

TEST(SparseMatrixTests, ResizeToSmaller1) {
  auto matrix = sparse<GMPRational>({
      {1, 2, 3},
      {4, 5, 6},
  });

  matrix.resize(1, 1);

  ASSERT_EQ(matrix.rows(), 1);
  ASSERT_EQ(matrix.cols(), 1);

  Matrix<GMPRational> expected = {{1}};

  ASSERT_EQ(Matrix(matrix), expected);
}

TEST(SparseMatrixTests, ResizeToSmaller2) {
  auto matrix = sparse<GMPRational>({
      {1, 2, 3},
      {4, 5, 6},
  });

  matrix.resize(1, 3);

  ASSERT_EQ(matrix.rows(), 1);
  ASSERT_EQ(matrix.cols(), 3);

  Matrix<GMPRational> expected = {{1, 2, 3}};

  ASSERT_EQ(Matrix(matrix), expected);
}

TEST(SparseMatrixTests, PrunesEntries1) {
  auto matrix = CSCMatrix<GMPRational>::zeros(10);

  matrix.add_column();
  matrix.push_to_last_column(1, 5);

  matrix.add_column();
  matrix.push_to_last_column(5, 5);

  matrix.resize(5, 5);

  ASSERT_EQ(matrix.rows(), 5);
  ASSERT_EQ(matrix.cols(), 5);

  ASSERT_EQ(matrix.entries_count(), 1);
}

TEST(SparseMatrixTests, PrunesEntries2) {
  auto matrix = CSCMatrix<GMPRational>::zeros(10);

  matrix.add_column();
  matrix.push_to_last_column(1, 5);

  matrix.add_column();
  matrix.push_to_last_column(5, 5);

  matrix.resize(20, 1);

  ASSERT_EQ(matrix.rows(), 20);
  ASSERT_EQ(matrix.cols(), 1);

  ASSERT_EQ(matrix.entries_count(), 1);
}
