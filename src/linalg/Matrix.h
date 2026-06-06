#pragma once

#include <format>
#include <vector>

#include "Concepts.h"
#include "expr/SubColsExpr.h"
#include "expr/SubRowsExpr.h"
#include "expr/TransposedExpr.h"

#include "Arithmetics.h"

namespace linalg {

template <typename Field>
class Matrix {
  size_t rows_;
  size_t cols_;
  std::vector<Field> data_;

  // protected so that Vector can access them
 protected:
  size_t get_index(size_t row, size_t col) const { return row * cols_ + col; }

  explicit Matrix(size_t rows, size_t cols)
      : rows_(rows), cols_(cols), data_(rows * cols) {}

  Matrix(size_t rows, size_t cols, Field value)
      : rows_(rows), cols_(cols), data_(rows * cols, value) {}

 public:
  using FieldType = Field;

  //
  Matrix() : Matrix(0, 0) {}

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

  template <MatrixRange T>
    requires std::same_as<MatrixFieldType<T>, Field>
  Matrix(T&& other) : Matrix(other.rows(), other.cols(), 0) {
    // TODO: check that i, j don't go outside of range
    for (const auto [i, j, value] : other.entries()) {
      (*this)[i, j] += value;
    }
  }

  Matrix(const Matrix& other) = default;

  Matrix(Matrix&& other)
      : rows_(other.rows_), cols_(other.cols_), data_(std::move(other.data_)) {
    other.rows_ = 0;
    other.cols_ = 0;
    other.data_.clear();
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

  template <MatrixRange T>
    requires std::same_as<MatrixFieldType<T>, Field>
  Matrix& operator=(T&& other) {
    // TODO: aliasing
    // TODO: check that i, j don't go outside of range
    rows_ = other.rows();
    cols_ = other.cols();

    data_.resize(rows_ * cols_);
    std::ranges::fill_n(data_, rows_ * cols_, 0);

    for (const auto [i, j, value] : other.entries()) {
      (*this)[i, j] += value;
    }

    return *this;
  }

  Matrix& operator=(const Matrix& other) {
    if (this != &other) {
      auto copy = other;
      std::swap(*this, copy);
    }

    return *this;
  }

  Matrix& operator=(Matrix&& other) {
    auto copy = std::move(other);
    std::swap(*this, copy);

    return *this;
  }

  //
  Field& operator[](size_t row, size_t col) {
    return data_[get_index(row, col)];
  }

  const Field& operator[](size_t row, size_t col) const {
    return data_[get_index(row, col)];
  }

  auto col_entries(size_t col) const {
    return std::views::iota(size_t{0}, rows()) |
           std::views::transform([this, col](size_t row) {
             return std::pair{row, (*this)[row, col]};
           });
  }

  auto row_entries(size_t row) const {
    return std::views::iota(size_t{0}, cols()) |
           std::views::transform([this, row](size_t col) {
             return std::pair{col, (*this)[row, col]};
           });
  }

  auto entries() {
    return std::views::iota(size_t{0}, data_.size()) |
           std::views::transform([this](size_t idx) {
             return std::tuple{idx / cols_, idx % cols_, data_[idx]};
           });
  }

  auto entries() const {
    return std::views::iota(size_t{0}, data_.size()) |
           std::views::transform([this](size_t idx) {
             return std::tuple{idx / cols_, idx % cols_, data_[idx]};
           });
  }

  template <IndicesRange RowRange, IndicesRange ColRange>
  auto operator[](RowRange&& rows, ColRange&& cols) const {
    return detail::SubRowsExpr(
        detail::SubColsExpr(*this, std::forward<ColRange>(cols)),
        std::forward<RowRange>(rows));
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

  auto transposed() const { return detail::TransposedExpr(*this); }

  //
  template <MatrixRange R>
  Matrix& operator+=(const R& other) {
    if (shape() != other.shape()) {
      throw std::invalid_argument(std::format(
          "Incompatible operand shape: {} != {}", shape(), other.shape()));
    }

    for (const auto [i, j, value] : other.entries()) {
      (*this)[i, j] += value;
    }

    return *this;
  }

  template <MatrixRange R>
  Matrix& operator-=(const R& other) {
    if (shape() != other.shape()) {
      throw std::invalid_argument(std::format(
          "Incompatible operand shape: {} != {}", shape(), other.shape()));
    }

    for (const auto [i, j, value] : other.entries()) {
      (*this)[i, j] -= value;
    }

    return *this;
  }
};

template <MatrixRange T>
Matrix(T&&) -> Matrix<MatrixFieldType<T>>;

}  // namespace linalg
