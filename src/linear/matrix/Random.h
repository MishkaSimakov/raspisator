#pragma once

#include <concepts>
#include <random>

#include "Matrix.h"
#include "Rank.h"

namespace linalg {

template <typename Field, std::invocable G>
Matrix<Field> random(size_t height, size_t width, G&& generator) {
  Matrix<Field> result(height, width, 0);

  for (size_t i = 0; i < height; ++i) {
    for (size_t j = 0; j < width; ++j) {
      result[i, j] = generator();
    }
  }

  return result;
}

template <typename Field, std::invocable G>
Matrix<Field> random_invertible(size_t size, G&& generator) {
  Matrix<Field> result(size, size, 0);

  do {
    result = random<Field>(size, size, std::forward<G>(generator));
  } while (linalg::rank(result) != size);

  return result;
}

}  // namespace linalg
