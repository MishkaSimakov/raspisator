#pragma once

#include <span>

#include "linalg/Matrix.h"
#include "linalg/Vector.h"
#include "utils/MaybeConst.h"

namespace linalg {

enum class EtaMatrixType { ROW, COLUMN };

// Apply methods can accept both Vector and Matrix. Vector has constexpr cols()
// method so that compiler can optimize out columns loop for Vector
// instantiation.
template <typename Field, bool is_const>
class SparseEtaMatrixView {
  using ValueT = maybe_const_t<std::pair<size_t, Field>, is_const>;

  size_t size_;
  size_t pivot_index_;
  EtaMatrixType type_;
  std::span<ValueT> entries_;

  // Returns element on diagonal in pivot row or column.
  Field pivot_diagonal() const {
    for (const auto& [row, value] : entries_) {
      if (row == pivot_index_) {
        return value;
      }
    }

    return 0;
  }

  template <typename T>
  void apply_as_col(T& matrix) const {
    for (size_t col = 0; col < matrix.cols(); ++col) {
      const Field a = matrix[pivot_index_, col];
      matrix[pivot_index_, col] = 0;

      for (const auto& [row, value] : entries_) {
        matrix[row, col] += a * value;
      }
    }
  }

  template <typename T>
  void apply_as_row(T& matrix) const {
    for (size_t col = 0; col < matrix.cols(); ++col) {
      Field dot = 0;

      for (const auto& [row, value] : entries_) {
        dot += value * matrix[row, col];
      }

      matrix[pivot_index_, col] = dot;
    }
  }

  template <typename T>
  void apply_inverse_as_col(T& matrix) const {
    for (size_t col = 0; col < matrix.cols(); ++col) {
      const Field a = matrix[pivot_index_, col];
      const Field diagonal = pivot_diagonal();

      for (const auto& [row, value] : entries_) {
        if (row == pivot_index_) {
          matrix[row, col] = a / diagonal;
        } else {
          matrix[row, col] -= a * value / diagonal;
        }
      }
    }
  }

  template <typename T>
  void apply_inverse_as_row(T& matrix) const {
    for (size_t col = 0; col < matrix.cols(); ++col) {
      Field dot = 0;

      for (auto [row, value] : entries_) {
        dot +=
            row == pivot_index_ ? matrix[row, col] : -value * matrix[row, col];
      }

      matrix[pivot_index_, col] = dot / pivot_diagonal();
    }
  }

 public:
  using FieldType = Field;

  SparseEtaMatrixView(size_t size, size_t pivot_index, EtaMatrixType type,
                      std::span<ValueT> values)
      : size_(size), pivot_index_(pivot_index), type_(type), entries_(values) {}

  // NOLINTNEXTLINE(google-explicit-constructor)
  SparseEtaMatrixView(SparseEtaMatrixView<Field, false> other)
    requires(is_const)
      : pivot_index_(other.pivot_index_),
        type_(other.type),
        entries_(other.entries_) {}

  template <typename T>
    requires std::same_as<T, Matrix<Field>> || std::same_as<T, Vector<Field>>
  T apply(T matrix) const {
    if (type_ == EtaMatrixType::ROW) {
      apply_as_row(matrix);
    } else {
      apply_as_col(matrix);
    }

    return std::move(matrix);
  }

  template <typename T>
    requires std::same_as<T, Matrix<Field>> || std::same_as<T, Vector<Field>>
  T apply_transposed(T matrix) const {
    if (type_ == EtaMatrixType::ROW) {
      apply_as_col(matrix);
    } else {
      apply_as_row(matrix);
    }

    return std::move(matrix);
  }

  template <typename T>
    requires std::same_as<T, Matrix<Field>> || std::same_as<T, Vector<Field>>
  T apply_inverse(T matrix) const {
    if (type_ == EtaMatrixType::ROW) {
      apply_inverse_as_row(matrix);
    } else {
      apply_inverse_as_col(matrix);
    }

    return std::move(matrix);
  }

  template <typename T>
    requires std::same_as<T, Matrix<Field>> || std::same_as<T, Vector<Field>>
  T apply_inverse_transposed(T matrix) const {
    if (type_ == EtaMatrixType::ROW) {
      apply_inverse_as_col(matrix);
    } else {
      apply_inverse_as_row(matrix);
    }

    return std::move(matrix);
  }

  //
  size_t pivot_index() const { return pivot_index_; }
  std::span<ValueT> pivot_entries() const { return entries_; }
  EtaMatrixType type() const { return type_; }

  //
  template <typename F>
  void entries(F&& f) const {
    for (size_t i = 0; i < size_; ++i) {
      if (i != pivot_index_) {
        f(i, i, Field(1));
      }
    }

    for (const auto& [i, value] : entries_) {
      if (type_ == EtaMatrixType::ROW) {
        f(pivot_index_, i, value);
      } else {
        f(i, pivot_index_, value);
      }
    }
  }

  //
  size_t rows() const { return size_; }
  size_t cols() const { return size_; }
  std::pair<size_t, size_t> shape() const { return {rows(), cols()}; }

  //
  Field det() const { return pivot_diagonal(); }
};

}  // namespace linalg
