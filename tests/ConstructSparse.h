#pragma once

#include "linalg/CSCMatrix.h"
#include "linalg/Matrix.h"

// Creates an integer CSCMatrix using dense initializer_list representation.
// Stores only non-zeros in the returned sparse matrix.
template <typename Field>
linalg::CSCMatrix<Field> sparse(
    std::initializer_list<std::initializer_list<Field>> il) {
  linalg::Matrix<Field> dense = il;

  auto sparse = linalg::CSCMatrix<Field>::zeros(dense.rows());

  for (size_t col = 0; col < dense.cols(); ++col) {
    sparse.add_column();

    for (size_t row = 0; row < dense.rows(); ++row) {
      if (dense[row, col] != 0) {
        sparse.push_to_last_column(row, dense[row, col]);
      }
    }
  }

  return sparse;
}
