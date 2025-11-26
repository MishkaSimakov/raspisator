#pragma once

#include <algorithm>

#include "linear/matrix/Matrix.h"
#include "linear/matrix/RowBasis.h"
#include "linear/model/LP.h"

template <typename Field>
bool check_bfs(const Matrix<Field>& A, const Matrix<Field>& b,
               const BFS<Field>& bfs) {
  if (A * bfs.point != b) {
    std::cout << A * bfs.point - b << std::endl;
    return false;
  }

  size_t n = A.get_height();

  if (bfs.basic_variables.size() != n) {
    return false;
  }

  Matrix<Field> A_b(n, n);

  for (size_t i = 0; i < n; ++i) {
    A_b[{0, n}, i] = A[{0, n}, bfs.basic_variables[i]];
  }

  // rk(A) == n
  auto row_basis = linalg::get_row_basis(A_b);
  if (row_basis.size() != n) {
    return false;
  }

  for (size_t i = 0; i < bfs.point.get_height(); ++i) {
    if (std::ranges::find(bfs.basic_variables, i) ==
        bfs.basic_variables.end()) {
      if (bfs.point[i, 0] != 0) {
        return false;
      }
    } else {
      if (bfs.point[i, 0] < 0) {
        return false;
      }
    }
  }

  return true;
}
