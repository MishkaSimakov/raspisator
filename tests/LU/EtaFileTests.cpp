#include <gtest/gtest.h>

#include "linear/BigInteger.h"
#include "linear/sparse/EtaFile.h"

TEST(EtaFileTests, ApplyInverseTest1) {
  linalg::EtaFile<Rational> file;

  const Matrix<Rational> values = {{1}, {2}, {3}};
  file.push_back(1, values, linalg::EtaType::COLUMN);

  const auto matrix = Matrix<Rational>::unity(3);
  const auto result = file.apply_inverse(matrix, *file.begin());

  const Matrix<Rational> expected = {
      {1, -Rational{1} / 2, 0},
      {0, Rational{1} / 2, 0},
      {0, -Rational{3} / 2, 1},
  };

  ASSERT_EQ(result, expected);
}

TEST(EtaFileTests, ApplyInverseTest2) {
  linalg::EtaFile<Rational> file;

  const Matrix<Rational> values = {{1}, {2}, {3}};
  file.push_back(1, values, linalg::EtaType::ROW);

  const auto matrix = Matrix<Rational>::unity(3);
  const auto result = file.apply_inverse(matrix, *file.begin());

  const Matrix<Rational> expected = {
      {1, 0, 0},
      {-Rational{1} / 2, Rational{1} / 2, -Rational{3} / 2},
      {0, 0, 1},
  };

  ASSERT_EQ(result, expected);
}
