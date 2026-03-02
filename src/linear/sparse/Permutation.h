#pragma once

#include <cassert>
#include <numeric>
#include <set>
#include <vector>

#include "linear/matrix/Matrix.h"
#include "linear/sparse/CSCMatrix.h"

namespace linalg {

// Stores permutation matrix P
class Permutation {
  // permutation_[i] = j means that i-th row in A is j-th row in PA
  std::vector<size_t> permutation_;

  explicit Permutation(std::vector<size_t> permutation)
      : permutation_(std::move(permutation)) {}

 public:
  static Permutation id(size_t size) {
    std::vector<size_t> permutation(size);
    std::iota(permutation.begin(), permutation.end(), 0);

    return Permutation(std::move(permutation));
  }

  static Permutation from_vector(std::vector<size_t> permutation) {
    // check permutation correctness (only in debug)
#ifndef NDEBUG
    std::set unique(permutation.begin(), permutation.end());
    assert(unique.size() == permutation.size());

    for (size_t element : permutation) {
      assert(element < permutation.size());
    }
#endif

    return Permutation(std::move(permutation));
  }

  // Returns PA - row permutation of A
  // Allocates new matrix of the same shape as A.
  // TODO: this may be implemented without extra allocations:
  // https://blog.merovius.de/posts/2014-08-12-applying-permutation-in-constant/
  template <typename Field>
  Matrix<Field> apply(const Matrix<Field>& A) const {
    auto [n, d] = A.shape();

    assert(n == size());

    auto result = Matrix<Field>(n, d);

    for (size_t row = 0; row < n; ++row) {
      for (size_t col = 0; col < d; ++col) {
        result[permutation_[row], col] = A[row, col];
      }
    }

    return result;
  }

  // Returns PA - row permutation of A
  // Does not allocate new matrices, performs transformation in place
  template <typename Field>
  CSCMatrix<Field> apply(CSCMatrix<Field> A) const {
    auto [n, d] = A.shape();

    assert(n == size());

    for (size_t col = 0; col < d; ++col) {
      for (size_t& row : A.get_column(col) | std::views::keys) {
        row = permutation_[row];
      }
    }

    return A;
  }

  // Returns PA - row permutation of A
  // Does not allocate new matrices, performs transformation in place
  template <typename Field>
  std::vector<std::pair<size_t, Field>> apply(
      std::vector<std::pair<size_t, Field>> A) const {
    assert(A.size() == size());

    for (size_t& row : A | std::views::keys) {
      row = permutation_[row];
    }

    return A;
  }

  // Returns P^T A - row permutation of A
  // Allocates new matrix of the same shape as A.
  template <typename Field>
  Matrix<Field> apply_transposed(const Matrix<Field>& A) const {
    auto [n, d] = A.shape();

    assert(n == size());

    auto result = Matrix<Field>(n, d);

    for (size_t row = 0; row < n; ++row) {
      for (size_t col = 0; col < d; ++col) {
        result[row, col] = A[permutation_[row], col];
      }
    }

    return result;
  }

  size_t size() const { return permutation_.size(); }

  // Returns the row index in PA where row `row` of A is mapped by P
  size_t operator[](size_t row) const {
    assert(row < size() && "row index out of bounds");

    return permutation_[row];
  }

  template <typename Field>
  Matrix<Field> as_dense_matrix() const {
    size_t n = size();
    Matrix<Field> result(n, n, 0);

    for (size_t row = 0; row < n; ++row) {
      result[permutation_[row], row] = 1;
    }

    return result;
  }

  template <typename Field>
  CSCMatrix<Field> as_sparse_matrix() const {
    size_t n = size();
    CSCMatrix<Field> result(n);

    for (size_t col = 0; col < n; ++col) {
      result.add_column();
      result.push_to_last_column(permutation_[col], Field(1));
    }

    return result;
  }

  Permutation transposed() const {
    size_t n = size();
    std::vector<size_t> result(n);

    for (size_t i = 0; i < n; ++i) {
      result[permutation_[i]] = i;
    }

    return from_vector(std::move(result));
  }

  // Swaps rows i and j of the permutation matrix P
  void swap(size_t i, size_t j) {
    assert(i < size() && "row index i out of bounds");
    assert(j < size() && "row index j out of bounds");

    std::swap(permutation_[i], permutation_[j]);
  }
};

}  // namespace linalg
