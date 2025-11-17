#pragma once

#include <iostream>
#include <sstream>
#include <vector>
#include <iomanip>

template <typename Field>
class Matrix {
  std::vector<Field> data_;
  size_t rows_count_;
  size_t cols_count_;

  size_t get_index(size_t row, size_t col) const {
    return col + cols_count_ * row;
  }

  void swap_rows(size_t first, size_t second) {
    if (first == second) {
      return;
    }

    for (size_t j = 0; j < cols_count_; ++j) {
      std::swap((*this)[first, j], (*this)[second, j]);
    }
  }

 public:
  Matrix(size_t rows, size_t cols)
      : data_(rows * cols), rows_count_(rows), cols_count_(cols) {}

  Matrix(size_t rows, size_t cols, const Field& value)
      : data_(rows * cols, value), rows_count_(rows), cols_count_(cols) {}

  Matrix(const std::initializer_list<std::initializer_list<Field>>& values)
      : Matrix(values.size(), values.begin()->size()) {
    size_t row = 0;
    size_t col = 0;

    for (const auto& row_value : values) {
      if (row_value.size() != cols_count_) {
        throw std::invalid_argument("All rows must have equal length.");
      }

      for (const auto& value : row_value) {
        data_[get_index(row, col)] = Field(value);
        ++col;
      }

      ++row;
      col = 0;
    }
  }

  Matrix get_extended(size_t new_height, size_t new_width,
                      Field fill_value) const {
    if (new_height < rows_count_ || new_width < cols_count_) {
      throw std::invalid_argument("New size must be larger than the old size.");
    }

    auto result = Matrix(new_height, new_width, fill_value);

    for (size_t i = 0; i < rows_count_; ++i) {
      for (size_t j = 0; j < cols_count_; ++j) {
        result[i, j] = (*this)[i, j];
      }
    }

    return result;
  }

  template <typename OtherField>
  static Matrix zeros_like(const Matrix<OtherField>& other) {
    return Matrix(other.shape().first, other.shape().second, 0);
  }

  template <typename Functor>
  static Matrix elementwise(const Matrix& left, const Matrix& right,
                            Functor functor) {
    if (left.shape() != right.shape()) {
      throw std::invalid_argument("Arguments shapes are not compatible.");
    }

    auto [n, m] = left.shape();
    auto result = Matrix(n, m);

    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < m; ++j) {
        result[i, j] = functor(left[i, j], right[i, j]);
      }
    }

    return result;
  }

  static Matrix unity(size_t size) {
    auto result = Matrix(size, size);

    for (size_t i = 0; i < size; ++i) {
      result[i, i] = 1;  // TODO: Field?
    }

    return result;
  }

  size_t get_width() const { return cols_count_; }
  size_t get_height() const { return rows_count_; }

  std::pair<size_t, size_t> shape() const { return {rows_count_, cols_count_}; }

  Field& operator[](size_t row, size_t col) {
    if (row >= rows_count_ || col >= cols_count_) {
      throw std::invalid_argument("The index is outside of matrix");
    }

    return data_[get_index(row, col)];
  }

  const Field& operator[](size_t row, size_t col) const {
    if (row >= rows_count_ || col >= cols_count_) {
      throw std::invalid_argument("The index is outside of matrix");
    }

    return data_[get_index(row, col)];
  }

  Matrix& inverse() {
    auto [n, m] = shape();
    if (n != m) {
      throw std::invalid_argument("Matrix must be square for inverse.");
    }

    auto inv = Matrix::unity(n);

    struct MultiplicationHook {
      Matrix& inv;

      explicit MultiplicationHook(Matrix& inv) : inv(inv) {}

      void operator()(size_t multiplicand, Field multiplier) {
        for (size_t j = 0; j < inv.shape().second; ++j) {
          inv[multiplicand, j] *= multiplier;
        }
      }
    };

    struct SubtractionHook {
      Matrix& inv;

      explicit SubtractionHook(Matrix& inv) : inv(inv) {}

      void operator()(size_t minuend, size_t subtrahend, Field multiplier) {
        for (size_t j = 0; j < inv.shape().second; ++j) {
          inv[minuend, j] -= inv[subtrahend, j] * multiplier;
        }
      }
    };

    for (size_t j = 0; j < n; ++j) {
      size_t nonzero_row_index = n;

      for (size_t i = j; i < n; ++i) {
        if ((*this)[i, j] != 0) {
          nonzero_row_index = i;
          break;
        }
      }

      if (nonzero_row_index == n) {
        throw std::invalid_argument("Matrix must be non-singular.");
      }

      inv.swap_rows(j, nonzero_row_index);
      swap_rows(j, nonzero_row_index);

      gaussian_elimination(j, j, SubtractionHook{inv}, MultiplicationHook{inv});
    }

    std::swap(inv, *this);

    return *this;
  }

  bool operator==(const Matrix& other) const {
    if (other.rows_count_ != rows_count_ || other.cols_count_ != cols_count_) {
      return false;
    }

    for (size_t i = 0; i < rows_count_; ++i) {
      for (size_t j = 0; j < cols_count_; ++j) {
        if ((*this)[i, j] != other[i, j]) {
          return false;
        }
      }
    }

    return true;
  }

  // gaussian elimination
  struct NoopSubtractHook {
    void operator()(size_t minuend, size_t subtrahend, Field multiplier) {}
  };

  struct NoopMultiplyHook {
    void operator()(size_t multiplicand, Field multiplier) {}
  };

  // subtracts row #row_index from other rows, so the column #col_index becomes
  // (0, ..., 1, 0, ..., 0)^T, where 1 is in row #row_index
  template <typename SubtractHook = NoopSubtractHook,
            typename MultiplyHook = NoopMultiplyHook>
  void gaussian_elimination(size_t row_index, size_t col_index,
                            SubtractHook subtract_hook = {},
                            MultiplyHook multiply_hook = {}) {
    if ((*this)[row_index, col_index] == 0) {
      throw std::invalid_argument(
          "Element must be non-zero for gaussian elimination.");
    }

    Field inverse = 1 / (*this)[row_index, col_index];
    multiply_hook(row_index, inverse);
    for (size_t j = 0; j < cols_count_; ++j) {
      (*this)[row_index, j] *= inverse;
    }

    for (size_t i = 0; i < rows_count_; ++i) {
      if (i == row_index || (*this)[i, col_index] == 0) {
        continue;
      }

      Field coef = (*this)[i, col_index];
      subtract_hook(i, row_index, coef);

      for (size_t j = 0; j < cols_count_; ++j) {
        (*this)[i, j] -= (*this)[row_index, j] * coef;
      }
    }
  }
};

template <typename Field>
Matrix<Field> operator*(const Matrix<Field>& left, const Matrix<Field>& right) {
  auto [n1, m1] = left.shape();
  auto [n2, m2] = right.shape();

  if (m1 != n2) {
    throw std::invalid_argument("Arguments shapes are not compatible.");
  }

  auto result = Matrix<Field>(n1, m2);

  for (size_t i = 0; i < n1; ++i) {
    for (size_t j = 0; j < m2; ++j) {
      for (size_t k = 0; k < n2; ++k) {
        result[i, j] += left[i, k] * right[k, j];
      }
    }
  }

  return result;
}

template <typename Field>
Matrix<Field> operator-(const Matrix<Field>& left, const Matrix<Field>& right) {
  return Matrix<Field>::elementwise(left, right,
                                    [](Field l, Field r) { return l - r; });
}

template <typename Field>
Matrix<Field> operator+(const Matrix<Field>& left, const Matrix<Field>& right) {
  return Matrix<Field>::elementwise(left, right,
                                    [](Field l, Field r) { return l + r; });
}

template <typename Field>
std::ostream& operator<<(std::ostream& os, const Matrix<Field>& matrix) {
  auto [n, m] = matrix.shape();

  std::vector<std::string> results(n * m);
  size_t max_length = 0;

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      std::stringstream ss;
      ss << matrix[i, j];
      results[i * m + j] = ss.str();

      max_length = std::max(max_length, results[i * m + j].size());
    }
  }

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      os << std::right << std::setw(max_length) << results[i * m + j] << " ";
    }

    os << "\n";
  }

  return os;
}
