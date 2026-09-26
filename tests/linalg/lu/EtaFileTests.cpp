#include <gtest/gtest.h>

#include "../../../src/field/BigInteger.h"
#include "linalg/Matrix.h"
#include "linalg/Print.h"
#include "linalg/lu/EtaFile.h"

using namespace linalg;

static_assert(MatrixRange<SparseEtaMatrixView<int, true>>);

TEST(EtaFileTests, ApplyInverseTest1) {
  EtaFile<Rational> file(3);

  const Vector<Rational> values = {1, 2, 3};
  file.push_back(1, values, EtaMatrixType::COLUMN);

  Vector<Rational> vector = {0, 1, 0};
  vector = (*file.begin()).apply_inverse(std::move(vector));

  Vector expected = {-Rational{1} / 2, Rational{1} / 2, -Rational{3} / 2};

  std::cout << *file.begin() << std::endl;

  ASSERT_EQ(vector, expected);
}

TEST(EtaFileTests, ApplyInverseTest2) {
  EtaFile<Rational> file(3);

  const Vector<Rational> values = {1, 2, 3};
  file.push_back(1, values, EtaMatrixType::ROW);

  Vector<Rational> vector = {1, 0, 0};
  vector = (*file.begin()).apply_inverse(std::move(vector));

  const Vector<Rational> expected = {1, -Rational{1} / 2, 0};

  ASSERT_EQ(vector, expected);
}
