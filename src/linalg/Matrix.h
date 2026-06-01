#pragma once

#include <format>
#include <vector>

#include "Types.h"

namespace linalg {

template <typename Field>
class Matrix {
  size_t rows_;
  size_t cols_;
  std::vector<Field> data_;

  size_t get_index(size_t row, size_t col) const { return row * cols_ + col; }

  explicit Matrix(size_t rows = 0, size_t cols = 0)
      : rows_(rows), cols_(cols), data_(rows * cols) {}

  Matrix(size_t rows, size_t cols, Field value)
      : rows_(rows), cols_(cols), data_(rows * cols, value) {}

 public:
  //
  Matrix(std::initializer_list<std::initializer_list<Field>> values)
      : Matrix(values.size(), values.begin()->size()) {
    size_t row = 0;
    size_t col = 0;

    for (const auto& row_value : values) {
      if (row_value.size() != cols_) {
        throw std::invalid_argument(
            std::format("Row {} has wrong size. Expected {}, but size is {}.",
                        row, cols_, row_value.size()));
      }

      for (const auto& value : row_value) {
        data_[get_index(row, col)] = value;
        ++col;
      }

      ++row;
      col = 0;
    }
  }

  static Matrix uninitialized(size_t rows, size_t cols) {
    return Matrix(rows, cols);
  }

  static Matrix zeros(size_t rows, size_t cols) {
    return Matrix(rows, cols, 0);
  }

  template <typename G>
    requires std::invocable<G, size_t, size_t> &&
             std::convertible_to<std::invoke_result_t<G, size_t, size_t>, Field>
  static Matrix generate(size_t rows, size_t cols, G&& generator) {
    auto result = Matrix(rows, cols);

    for (size_t i = 0; i < rows; ++i) {
      for (size_t j = 0; j < cols; ++j) {
        result[i, j] = generator(i, j);
      }
    }

    return result;
  }

  //
  Field& operator[](size_t row, size_t col) {
    return data_[get_index(row, col)];
  }

  const Field& operator[](size_t row, size_t col) const {
    return data_[get_index(row, col)];
  }

  template <IndicesRange RowRange, IndicesRange ColRange>
  auto operator[](const RowRange& rows, const ColRange& cols) {
    return IndexedView<Field, Matrix, RowRange, ColRange>{*this, rows, cols};
  }

  template <IndicesRange RowRange, IndicesRange ColRange>
  auto operator[](const RowRange& rows, const ColRange& cols) const {
    return IndexedView<Field, Matrix, RowRange, ColRange>{*this, rows, cols};
  }

  //
  size_t rows() const { return rows_; }
  size_t cols() const { return cols_; }
  std::pair<size_t, size_t> shape() const { return {rows(), cols()}; }

  //
  void resize(size_t new_rows, size_t new_cols) {
    if (new_rows >= rows_ && new_cols >= cols_) {
      data_.resize(new_rows * new_cols, 0);

      for (size_t i = 0; i < rows_; ++i) {
        for (size_t j = 0; j < new_cols; ++j) {
          const size_t row = rows_ - i - 1;
          const size_t col = new_cols - j - 1;

          data_[row * new_cols + col] =
              col < cols_ ? data_[row * cols_ + col] : 0;
        }
      }

      rows_ = new_rows;
      cols_ = new_cols;
    } else if (new_rows >= rows_) {
      for (size_t row = 0; row < rows_; ++row) {
        for (size_t col = 0; col < new_cols; ++col) {
          data_[row * new_cols + col] = data_[row * cols_ + col];
        }
      }

      data_.resize(new_rows * new_cols);

      for (size_t row = rows_; row < new_rows; ++row) {
        for (size_t col = 0; col < new_cols; ++col) {
          data_[row * new_cols + col] = 0;
        }
      }

      rows_ = new_rows;
      cols_ = new_cols;
    } else if (new_cols >= cols_) {
      data_.resize(new_rows * new_cols, 0);

      for (size_t i = 0; i < new_rows; ++i) {
        for (size_t j = 0; j < new_cols; ++j) {
          const size_t row = new_rows - i - 1;
          const size_t col = new_cols - j - 1;

          data_[row * new_cols + col] =
              col < cols_ ? data_[row * cols_ + col] : 0;
        }
      }

      rows_ = new_rows;
      cols_ = new_cols;
    } else {
      for (size_t row = 0; row < new_rows; ++row) {
        for (size_t col = 0; col < new_cols; ++col) {
          data_[row * new_cols + col] = data_[row * cols_ + col];
        }
      }

      data_.resize(new_rows * new_cols);
      rows_ = new_rows;
      cols_ = new_cols;
    }
  }

  //
  bool operator==(const Matrix&) const = default;
};

}  // namespace linalg
