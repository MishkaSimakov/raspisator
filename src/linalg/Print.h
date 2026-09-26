#pragma once

#include <iomanip>
#include <iostream>

#include "Concepts.h"
#include "utils/Accumulators.h"

namespace linalg {

template <MatrixRange T>
std::ostream& operator<<(std::ostream& os, T&& matrix) {
  using std::to_string;
  using Field = MatrixFieldType<T>;

  auto [n, m] = matrix.shape();

  Matrix<Field> dense = matrix;
  std::vector<std::string> result(n * m, "-");

  Maximum<size_t> max_length;
  max_length.record(1);

  for (size_t row = 0; row < n; ++row) {
    for (size_t col = 0; col < m; ++col) {
      if constexpr (std::is_convertible_v<Field, std::string>) {
        result[row * m + col] = dense[row, col];
      } else {
        result[row * m + col] = to_string(dense[row, col]);
      }

      max_length.record(result[row * m + col].size());
    }
  }

  for (size_t i = 0; i < n; ++i) {
    os << "{";
    for (size_t j = 0; j < m; ++j) {
      os << std::right << std::setw(*max_length) << result[i * m + j];

      if (j + 1 != m) {
        os << ", ";
      }
    }

    os << "},\n";
  }

  return os;
}

}  // namespace linalg
