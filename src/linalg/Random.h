#pragma once

#include <concepts>
#include <random>

#include "CSCMatrix.h"
#include "Matrix.h"
#include "Rank.h"

#include "linalg/Print.h"

namespace linalg::random {

template <typename Field, typename G, typename D>
  requires std::invocable<D&, G&> &&
           std::uniform_random_bit_generator<std::remove_cvref_t<G>> &&
           std::convertible_to<std::invoke_result_t<D&, G&>, Field>
Matrix<Field> dense(size_t rows, size_t cols, G&& generator, D&& distribution) {
  auto result = Matrix<Field>::uninitialized(rows, cols);

  for (size_t row = 0; row < rows; ++row) {
    for (size_t col = 0; col < cols; ++col) {
      result[row, col] = distribution(generator);
    }
  }

  return result;
}

template <typename G, typename D>
auto dense(size_t rows, size_t cols, G&& generator, D&& distribution) {
  return dense<std::invoke_result_t<D&, G&>>(
      rows, cols, std::forward<G>(generator), std::forward<D>(distribution));
}

template <typename Field, typename G, typename D>
  requires std::invocable<D&, G&> &&
           std::uniform_random_bit_generator<std::remove_cvref_t<G>> &&
           std::convertible_to<std::invoke_result_t<D&, G&>, Field>
CSCMatrix<Field> sparse(size_t rows, size_t cols, size_t average_per_col,
                        G&& generator, D&& distribution) {
  auto result = CSCMatrix<double>::zeros(rows);

  std::bernoulli_distribution nonzero_distribution(
      static_cast<double>(average_per_col) / static_cast<double>(rows));

  for (size_t col = 0; col < cols; ++col) {
    result.add_column();

    for (size_t row = 0; row < rows; ++row) {
      if (nonzero_distribution(generator)) {
        result.push_to_last_column(row, distribution(generator));
      }
    }
  }

  return result;
}

template <typename G, typename D>
auto sparse(size_t rows, size_t cols, size_t average_per_col, G&& generator,
            D&& distribution) {
  return sparse<std::invoke_result_t<D&, G&>>(rows, cols, average_per_col,
                                              std::forward<G>(generator),
                                              std::forward<D>(distribution));
}

// Note: this method may be very inefficient for finite fields where the
// probability of singular matrix is considerable.
template <typename Field, typename G, typename D>
  requires std::invocable<D&, G&> &&
           std::uniform_random_bit_generator<std::remove_cvref_t<G>> &&
           std::convertible_to<std::invoke_result_t<D&, G&>, Field>
Matrix<Field> dense_invertible(size_t size, G&& generator, D&& distribution) {
  Matrix<Field> result;

  do {
    result = dense<Field>(size, size, generator, distribution);
  } while (rank(result) != size);

  return result;
}

template <typename G, typename D>
auto dense_invertible(size_t size, G&& generator, D&& distribution) {
  return dense_invertible<std::invoke_result_t<D&, G&>>(
      size, std::forward<G>(generator), std::forward<D>(distribution));
}

}  // namespace linalg::random
