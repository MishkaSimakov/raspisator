#pragma once

#include <cassert>
#include <numeric>
#include <set>
#include <utility>
#include <vector>

#include "CSCMatrix.h"
#include "Matrix.h"
#include "Vector.h"

namespace linalg {

// Stores permutation matrix P
class Permutation {
  static_assert(sizeof(size_t) == 8, "Different size_t size is not supported.");

  // flags
  constexpr static size_t loop_bit = static_cast<size_t>(1) << 63;
  constexpr static size_t visited_bit = static_cast<size_t>(1) << 62;

  constexpr static size_t flags = loop_bit | visited_bit;

  // permutation_[i] & ~flags = j means that i-th row in A is j-th row in PA.
  // Permutation can be decomposed into loops. (permutation_[i] & loop_bit) == 1
  // iff i is the smallest index in its loop.
  std::vector<size_t> permutation_;

  explicit Permutation(std::vector<size_t> permutation)
      : permutation_(std::move(permutation)) {
    for (size_t i = 0; i < permutation_.size(); ++i) {
      if ((permutation_[i] & visited_bit) != 0) {
        continue;
      }

      permutation_[i] = permutation_[i] | loop_bit;

      const size_t start = i;
      size_t current = i;

      do {
        permutation_[current] = permutation_[current] | visited_bit;

        current = permutation_[current] & ~flags;
      } while (current != start);
    }
  }

  template <typename T>
  struct MatrixOrVector : std::bool_constant<false> {};

  template <typename Field>
  struct MatrixOrVector<Vector<Field>> : std::bool_constant<true> {};

  template <typename Field>
  struct MatrixOrVector<Matrix<Field>> : std::bool_constant<true> {};

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
  template <typename T>
    requires MatrixOrVector<T>::value
  T apply(T A) const {
    auto [n, d] = A.shape();

    assert(n == size());

    for (size_t col = 0; col < d; ++col) {
      for (size_t row = 0; row < n; ++row) {
        if ((permutation_[row] & loop_bit) == 0) {
          continue;
        }

        const size_t start = row;
        size_t current = row;

        auto prev = A[row, col];

        do {
          auto next = A[(*this)[current], col];
          A[(*this)[current], col] = prev;
          prev = next;

          current = permutation_[current] & ~flags;
        } while (current != start);
      }
    }

    return std::move(A);
  }

  // Returns PA - row permutation of A
  template <typename Field>
  CSCMatrix<Field> apply(CSCMatrix<Field> A) const {
    auto [n, d] = A.shape();

    assert(n == size());

    for (size_t col = 0; col < d; ++col) {
      for (auto& [row, _] : A.get_column(col)) {
        row = (*this)[row];
      }
    }

    return A;
  }

  // Returns PA - row permutation of A
  template <typename Field>
  std::vector<std::pair<size_t, Field>> apply(
      std::vector<std::pair<size_t, Field>> A) const {
    assert(A.size() == size());

    for (size_t& row : A | std::views::keys) {
      row = (*this)[row];
    }

    return A;
  }

  // Returns P^T A - row permutation of A
  template <typename T>
    requires MatrixOrVector<T>::value
  T apply_transposed(T A) const {
    auto [n, d] = A.shape();

    assert(n == size());

    for (size_t col = 0; col < d; ++col) {
      for (size_t row = 0; row < n; ++row) {
        if ((permutation_[row] & loop_bit) == 0) {
          continue;
        }

        const size_t start = row;
        size_t current = row;

        auto old_value = A[row, col];

        while (true) {
          const size_t next = permutation_[current] & ~flags;

          if (next == start) {
            A[current, col] = old_value;
            break;
          }

          A[current, col] = A[next, col];
          current = next;
        }
      }
    }

    return std::move(A);
  }

  // Returns AP - column permutation of A
  // Allocates new matrix of the same shape as A.
  template <typename T>
    requires MatrixOrVector<T>::value
  T post_apply(T A) const {
    auto [n, d] = A.shape();

    assert(d == size());

    for (size_t row = 0; row < n; ++row) {
      for (size_t col = 0; col < d; ++col) {
        if ((permutation_[col] & loop_bit) == 0) {
          continue;
        }

        const size_t start = col;
        size_t current = col;

        auto old_value = A[row, col];

        while (true) {
          const size_t next = permutation_[current] & ~flags;

          if (next == start) {
            A[row, current] = old_value;
            break;
          }

          A[row, current] = A[row, next];
          current = next;
        }
      }
    }

    return std::move(A);
  }

  size_t size() const { return permutation_.size(); }

  // Returns the row index in PA where row `row` of A is mapped by P
  size_t apply(size_t row) const {
    assert(row < size() && "row index out of bounds");

    return permutation_[row] & ~flags;
  }

  size_t operator[](size_t row) const { return apply(row); }

  // Returns the column index in AP where column `col` of A is mapped by P
  size_t post_apply(size_t col) const {
    assert(col < size() && "col index out of bounds");

    size_t current = col;

    while (true) {
      const size_t next = permutation_[current] & ~flags;

      if (next == col) {
        return current;
      }

      current = next;
    }
  }

  template <typename Field>
  Matrix<Field> as_dense_matrix() const {
    size_t n = size();
    Matrix<Field> result(n, n, 0);

    for (size_t row = 0; row < n; ++row) {
      result[(*this)[row], row] = 1;
    }

    return result;
  }

  template <typename Field>
  CSCMatrix<Field> as_sparse_matrix() const {
    size_t n = size();
    auto result = CSCMatrix<Field>::zeros(n);

    for (size_t col = 0; col < n; ++col) {
      result.add_column();
      result.push_to_last_column((*this)[col], Field(1));
    }

    return result;
  }

  Permutation transposed() const {
    size_t n = size();
    std::vector<size_t> result(n);

    for (size_t i = 0; i < n; ++i) {
      result[(*this)[i]] = i;
    }

    return from_vector(std::move(result));
  }

  bool is_even() const {
    bool result = true;

    for (size_t i = 0; i < permutation_.size(); ++i) {
      if ((permutation_[i] & loop_bit) == 0) {
        continue;
      }

      const size_t start = i;
      size_t current = i;

      size_t size = 0;

      do {
        ++size;
        current = permutation_[current] & ~flags;
      } while (current != start);

      if (size % 2 == 0) {
        result = !result;
      }
    }

    return result;
  }
};

template <typename Field>
Matrix<Field> operator*(const Permutation& P, const Matrix<Field>& A) {
  return P.apply(A);
}

template <typename Field>
Matrix<Field> operator*(const Matrix<Field>& A, const Permutation& P) {
  return P.post_apply(A);
}

}  // namespace linalg
