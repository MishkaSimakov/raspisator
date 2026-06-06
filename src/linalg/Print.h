#pragma once

#include <iomanip>
#include <iostream>

#include "Concepts.h"
#include "utils/Accumulators.h"

namespace linalg {

template <SomeMatrixLike T>
std::ostream& operator<<(std::ostream& os, const T& matrix) {
  using std::to_string;

  auto [n, m] = matrix.shape();

  std::vector<std::string> result(n * m, "-");

  Maximum<size_t> max_length;
  max_length.record(1);

  for (const auto [row, col, value] : matrix.entries()) {
    if constexpr (std::is_convertible_v<typename T::FieldType, std::string>) {
      result[row * m + col] = value;
    } else {
      result[row * m + col] = to_string(value);
    }

    max_length.record(result[row * m + col].size());
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
