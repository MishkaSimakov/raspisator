#include <gtest/gtest.h>

#include "linalg/Matrix.h"
#include "linalg/Print.h"
#include "linalg/lu/EtaFile.h"
#include "support/GMPRational.h"

using namespace linalg;

static_assert(MatrixRange<SparseEtaMatrixView<int, true>>);

TEST(EtaFileTests, ApplyInverseTest1) {
  EtaFile<GMPRational> file(3);

  const Vector<GMPRational> values = {1, 2, 3};
  file.push_back(1, values, EtaMatrixType::COLUMN);

  Vector<GMPRational> vector = {0, 1, 0};
  vector = (*file.begin()).apply_inverse(std::move(vector));

  Vector expected = {-GMPRational{1} / 2, GMPRational{1} / 2,
                     -GMPRational{3} / 2};

  std::cout << *file.begin() << std::endl;

  ASSERT_EQ(vector, expected);
}

TEST(EtaFileTests, ApplyInverseTest2) {
  EtaFile<GMPRational> file(3);

  const Vector<GMPRational> values = {1, 2, 3};
  file.push_back(1, values, EtaMatrixType::ROW);

  Vector<GMPRational> vector = {1, 0, 0};
  vector = (*file.begin()).apply_inverse(std::move(vector));

  const Vector<GMPRational> expected = {1, -GMPRational{1} / 2, 0};

  ASSERT_EQ(vector, expected);
}
