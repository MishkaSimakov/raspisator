#pragma once

#include <iostream>
#include <ranges>
#include <span>
#include <vector>

#include "linear/FieldTraits.h"
#include "linear/matrix/Matrix.h"

// Matrix came to another matrix and said:
// - You know you are fat, aren't you?
// - What?
// - You are squarely filled with values.
// - What?
// - Oh, don't be so dense!

static size_t misses = 0;

template <typename Field>
class CSCMatrix {
  std::vector<Field> data_;
  std::vector<size_t> indices_;
  std::vector<size_t> index_pointers_;

  size_t rows_cnt_;

 public:
  // creates (height, 0) sparse matrix
  explicit CSCMatrix(size_t height)
      : index_pointers_(1, 0), rows_cnt_(height) {}

  explicit CSCMatrix(const Matrix<Field>& matrix)
      : index_pointers_(matrix.get_width() + 1),
        rows_cnt_(matrix.get_height()) {
    auto [n, d] = matrix.shape();

    index_pointers_[0] = 0;
    size_t nonzero_cnt = 0;

    for (size_t col = 0; col < d; ++col) {
      for (size_t row = 0; row < n; ++row) {
        if (FieldTraits<Field>::is_nonzero(matrix[row, col])) {
          data_.push_back(matrix[row, col]);
          indices_.push_back(row);
          ++nonzero_cnt;
        }
      }

      index_pointers_[col + 1] = nonzero_cnt;
    }
  }

  static CSCMatrix from_columns(const Matrix<Field>& dense,
                                const std::vector<size_t>& columns) {
    auto [n, d] = dense.shape();
    CSCMatrix result(n);

    for (size_t col : columns) {
      result.add_column();

      for (size_t row = 0; row < n; ++row) {
        if (FieldTraits<Field>::is_nonzero(dense[row, col])) {
          result.push_to_last_column(row, dense[row, col]);
        }
      }
    }

    return result;
  }

  std::pair<size_t, size_t> shape() const {
    return {rows_cnt_, index_pointers_.size() - 1};
  }

  auto get_column(size_t col) const {
    ssize_t begin = index_pointers_[col];
    ssize_t end = index_pointers_[col + 1];

    std::span<const size_t> indices(indices_.begin() + begin,
                                    indices_.begin() + end);
    std::span<const Field> values(data_.begin() + begin, data_.begin() + end);

    return std::views::zip(indices, values);
  }

  auto get_column(size_t col) {
    ssize_t begin = index_pointers_[col];
    ssize_t end = index_pointers_[col + 1];

    std::span<size_t> indices(indices_.begin() + begin, indices_.begin() + end);
    std::span<Field> values(data_.begin() + begin, data_.begin() + end);

    return std::views::zip(indices, values);
  }

  void add_column(std::span<const Field> dense, size_t shift = 0) {
    size_t height = rows_cnt_;

    if (dense.size() + shift > height) {
      throw std::invalid_argument("Column is too large for the matrix.");
    }

    size_t nonzero_cnt = index_pointers_.back();
    for (size_t i = 0; i < dense.size(); ++i) {
      if (FieldTraits<Field>::is_nonzero(dense[i])) {
        data_.push_back(dense[i]);
        indices_.push_back(i + shift);

        ++nonzero_cnt;
      }
    }

    index_pointers_.push_back(nonzero_cnt);
  }

  // adds zero column
  void add_column() { index_pointers_.push_back(index_pointers_.back()); }

  void push_to_last_column(size_t row, const Field& value) {
    if (!FieldTraits<Field>::is_nonzero(value)) {
      ++misses;
      return;
    }

    ++index_pointers_.back();

    indices_.push_back(row);
    data_.push_back(value);
  }

  // temporary
  void swap_rows(size_t i, size_t j) {
    for (size_t col = 0; col < (*this).shape().second; ++col) {
      for (size_t& row : get_column(col) | std::views::keys) {
        if (row == i) {
          row = j;
        } else if (row == j) {
          row = i;
        }
      }
    }
  }
};

template <typename Field>
std::ostream& operator<<(std::ostream& os, const CSCMatrix<Field>& matrix) {
  auto [n, m] = matrix.shape();

  Matrix<std::string> result(n, m, "-");

  for (size_t col = 0; col < m; ++col) {
    for (const auto& [row, value] : matrix.get_column(col)) {
      std::stringstream ss;
      ss << value;
      result[row, col] = ss.str();
    }
  }

  os << result;

  return os;
}

namespace linalg {

template <typename Field>
Matrix<Field> to_dense(const CSCMatrix<Field>& sparse) {
  auto [n, m] = sparse.shape();
  Matrix<Field> result(n, m);

  for (size_t col = 0; col < m; ++col) {
    for (const auto& [row, value] : sparse.get_column(col)) {
      result[row, col] = value;
    }
  }

  return result;
}

}  // namespace linalg
