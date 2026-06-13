#include <gtest/gtest.h>

#include "ConstructSparse.h"
#include "field/BigInteger.h"
#include "field/FieldTraits.h"
#include "linalg/Det.h"
#include "linalg/Linalg.h"
#include "linalg/Random.h"
#include "linalg/Stack.h"
#include "linalg/lu/LUPA.h"

using namespace linalg;

TEST(LUPATests, SmallSolveLinearTransposed) {
  const auto A = sparse<Rational>({
      {3, -7, -2, 2},
      {-3, 5, 1, 0},
      {6, -4, 0, -5},
      {-9, 5, -5, 12},
  });

  auto lupa = linalg::LUPA(A);

  lupa.set_columns(std::vector<size_t>{0, 1, 2, 3});

  Matrix<Rational> b = {{0}, {2}, {0}, {1}};
  auto solution = lupa.solve_linear_transposed(b);

  Matrix<Rational> expected = {{30}, {50}, {7}, {-2}};

  ASSERT_EQ(solution, expected);
}

TEST(LUPATests, SmallGetRow) {
  const auto A = sparse<Rational>({
      {3, -7, -2, 2},
      {-3, 5, 1, 0},
      {6, -4, 0, -5},
      {-9, 5, -5, 12},
  });

  auto lupa = linalg::LUPA(A);

  lupa.set_columns(std::vector<size_t>{0, 1, 2, 3});

  auto row = lupa.get_row(2);
  Matrix<Rational> expected = {{11}, {18}, {2}, {-1}};

  ASSERT_EQ(row, expected);
}

TEST(LUPATests, ChangeColumn) {
  const auto A = sparse<Rational>({
      {3, -7, -2, 2, 1, 1},
      {-3, 5, 1, 0, 0, 2},
      {6, -4, 0, -5, 2, 3},
      {-9, 5, -5, 12, 3, 4},
  });

  auto lupa = linalg::LUPA(A);

  lupa.set_columns(std::vector<size_t>{0, 1, 2, 3});

  lupa.change_column(1, 4);
  lupa.change_column(2, 5);

  auto inverse = lupa.get_inverse();

  const auto expected_cols = std::vector<size_t>{0, 4, 5, 3};
  auto expected = A.select_columns(expected_cols);

  ASSERT_EQ(inverse * expected, Matrix<Rational>::identity(4));
}

TEST(LUPATests, GetInverseMatrix) {
  const auto A = sparse<Rational>({
      {1, 0, 0},
      {0, 2, 1},
      {0, 1, 0},
  });

  auto lupa = linalg::LUPA(A);
  lupa.set_columns({0, 1, 2});

  const auto inverse = lupa.get_inverse();
  const Matrix<Rational> expected = {
      {1, 0, 0},
      {0, 0, 1},
      {0, 1, -2},
  };

  ASSERT_EQ(inverse, expected);
}

TEST(LUPATests, GetMatrix) {
  const auto A = sparse<Rational>({
      {1, 1, 0},
      {0, 2, 1},
      {0, 3, 0},
  });

  auto lupa = linalg::LUPA(A);
  lupa.set_columns({0, 1, 2});

  const auto matrix = lupa.get_matrix();

  ASSERT_EQ(matrix, A);
}

TEST(LUPATests, GetMatrixRandom) {
  constexpr size_t size = 10;

  std::default_random_engine random;
  std::uniform_int_distribution<int> value_distribution(-10, 10);

  for (size_t i = 0; i < 100; ++i) {
    auto A =
        random::dense_invertible<Rational>(size, random, value_distribution);
    auto sparse = CSCMatrix(A);

    auto lupa = linalg::LUPA(sparse);

    std::vector<size_t> columns(size);
    std::iota(columns.begin(), columns.end(), 0);
    lupa.set_columns(columns);

    const auto matrix = lupa.get_matrix();

    ASSERT_EQ(matrix, A);
  }
}

TEST(LUPATests, ChangeColumnsRandom) {
  constexpr size_t size = 10;

  std::default_random_engine random;
  std::uniform_int_distribution<int> value_distribution(-10, 10);

  for (size_t i = 0; i < 100; ++i) {
    SCOPED_TRACE(std::format("iteration: {}", i));

    const auto core =
        random::dense_invertible<Rational>(size, random, value_distribution);

    const auto dense = hstack(core, core);
    const auto sparse = CSCMatrix(dense);

    auto lupa = linalg::LUPA(sparse);

    std::vector<size_t> columns(size);
    std::iota(columns.begin(), columns.end(), 0);
    lupa.set_columns(columns);

    for (size_t j = 0; j < size; ++j) {
      lupa.change_column(j, j + size);
    }

    const auto matrix = lupa.get_matrix();

    ASSERT_EQ(matrix, core);
  }
}

TEST(LUPATests, ChangeColumnsRandomRoundtrip) {
  constexpr size_t size = 10;

  std::default_random_engine random;
  std::uniform_int_distribution<int> value_distribution(-10, 10);

  for (size_t i = 0; i < 100; ++i) {
    SCOPED_TRACE(std::format("iteration: {}", i));

    const auto core =
        random::dense_invertible<Rational>(size, random, value_distribution);
    const auto dense = hstack(core, core);
    const auto sparse = CSCMatrix(dense);

    auto lupa = linalg::LUPA(sparse);

    std::vector<size_t> columns(size);
    std::iota(columns.begin(), columns.end(), 0);
    lupa.set_columns(columns);

    for (size_t j = 0; j < size; ++j) {
      lupa.change_column(j, j + size);
    }

    for (size_t j = 0; j < size; ++j) {
      lupa.change_column(j, j);
    }

    const auto matrix = lupa.get_matrix();

    ASSERT_EQ(matrix, core);
  }
}

TEST(LUPATests, ChangeColumnsAndPurgeRandom) {
  constexpr size_t size = 10;

  std::default_random_engine random;
  std::uniform_int_distribution<int> value_distribution(-10, 10);

  for (size_t i = 0; i < 100; ++i) {
    SCOPED_TRACE(std::format("iteration: {}", i));

    const auto core =
        random::dense_invertible<Rational>(size, random, value_distribution);
    const auto dense = hstack(core, core);
    const auto sparse = CSCMatrix(dense);

    linalg::LUPAConfig config{
        .purge_after_iterations = 5,
        .refactorize_after_iterations = 10000,
    };

    auto lupa = linalg::LUPA(sparse, config);

    std::vector<size_t> columns(size);
    std::iota(columns.begin(), columns.end(), 0);
    lupa.set_columns(columns);

    for (size_t j = 0; j < size; ++j) {
      lupa.change_column(j, j + size);
    }

    for (size_t j = 0; j < size; ++j) {
      lupa.change_column(j, j);
    }

    const auto matrix = lupa.get_matrix();

    ASSERT_EQ(matrix, core);
  }
}

TEST(LUPATests, AccessWithoutSetColumnsThrowsOnFreshLupa) {
  const auto A = sparse<Rational>({
      {1, 2},
      {3, 4},
  });

  auto lupa = linalg::LUPA(A);  // set_columns never called

  const Vector<Rational> b = {1, 0};

  ASSERT_ANY_THROW(lupa.solve_linear(b));
  ASSERT_ANY_THROW(lupa.solve_linear_transposed(b));
  ASSERT_ANY_THROW(lupa.get_row(0));
  ASSERT_ANY_THROW(lupa.get_matrix());
  ASSERT_ANY_THROW(lupa.get_inverse());
  ASSERT_ANY_THROW(lupa.det());
}

// B = I_3, det(B) = 1
TEST(LUPATests, Det_IdentityMatrix) {
  const auto A = sparse<Rational>({
      {1, 0, 0},
      {0, 1, 0},
      {0, 0, 1},
  });

  auto lupa = linalg::LUPA(A);
  lupa.set_columns({0, 1, 2});

  const Rational expected = 1;

  ASSERT_EQ(lupa.det(), expected);
}

TEST(LUPATests, Det_SmallPositive) {
  const auto A = sparse<Rational>({
      {2, 1},
      {3, 4},
  });

  auto lupa = linalg::LUPA(A);
  lupa.set_columns({0, 1});

  ASSERT_EQ(lupa.det(), Rational{1} / 5);
}

TEST(LUPATests, Det_NegativeDeterminant) {
  const auto A = sparse<Rational>({
      {1, 2},
      {3, 4},
  });

  auto lupa = linalg::LUPA(A);
  lupa.set_columns({0, 1});

  ASSERT_EQ(lupa.det(), -Rational{1} / 2);
}

// Verify det matches naive_det when the initial basis is a non-trivial
// submatrix.
TEST(LUPATests, Det_SubmatrixColumnSelection) {
  const auto A = sparse<Rational>({
      {3, -7, -2, 2, 1, 1},
      {-3, 5, 1, 0, 0, 2},
      {6, -4, 0, -5, 2, 3},
      {-9, 5, -5, 12, 3, 4},
  });

  auto lupa = linalg::LUPA(A);
  lupa.set_columns({0, 4, 5, 3});

  const Matrix<Rational> B(A.select_columns(std::vector<size_t>{0, 4, 5, 3}));

  ASSERT_EQ(lupa.det(), 1 / det(B));
}

// det is correctly updated after a single column replacement.
TEST(LUPATests, Det_AfterChangeColumn) {
  const auto A = sparse<Rational>({
      {3, -7, -2, 2, 1, 1},
      {-3, 5, 1, 0, 0, 2},
      {6, -4, 0, -5, 2, 3},
      {-9, 5, -5, 12, 3, 4},
  });

  auto lupa = linalg::LUPA(A);
  lupa.set_columns({0, 1, 2, 3});

  lupa.change_column(1, 4);

  const Matrix<Rational> B(A.select_columns(std::vector<size_t>{0, 4, 2, 3}));

  ASSERT_EQ(lupa.det(), 1 / det(B));
}

// det remains correct after two sequential column replacements.
TEST(LUPATests, Det_AfterMultipleColumnChanges) {
  const auto A = sparse<Rational>({
      {3, -7, -2, 2, 1, 1},
      {-3, 5, 1, 0, 0, 2},
      {6, -4, 0, -5, 2, 3},
      {-9, 5, -5, 12, 3, 4},
  });

  auto lupa = linalg::LUPA(A);
  lupa.set_columns({0, 1, 2, 3});

  lupa.change_column(1, 4);
  lupa.change_column(2, 5);

  const Matrix<Rational> B(A.select_columns(std::vector<size_t>{0, 4, 5, 3}));

  ASSERT_EQ(lupa.det(), 1 / det(B));
}

// For random bases, det() matches naive_det before and after each column
// change.
TEST(LUPATests, Det_RandomBasisChanges) {
  constexpr size_t size = 5;

  std::default_random_engine random;
  std::uniform_int_distribution<int> value_distribution(-5, 5);

  for (size_t i = 0; i < 100; ++i) {
    SCOPED_TRACE(std::format("iteration: {}", i));

    const auto left =
        random::dense_invertible<Rational>(size, random, value_distribution);
    const auto dense = hstack(left, left);
    const auto sparse_A = CSCMatrix(dense);

    auto lupa = linalg::LUPA(sparse_A);
    std::vector<size_t> columns(size);
    std::iota(columns.begin(), columns.end(), 0);
    lupa.set_columns(columns);

    ASSERT_EQ(lupa.det(), 1 / det(left));

    for (size_t j = 0; j < size; ++j) {
      lupa.change_column(j, j + size);
      columns[j] = j + size;

      const Matrix<Rational> B(sparse_A.select_columns(columns));
      ASSERT_EQ(lupa.det(), 1 / det(B));
    }
  }
}
