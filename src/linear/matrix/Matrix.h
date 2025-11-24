#pragma once

#include <iomanip>
#include <span>
#include <sstream>
#include <vector>

#include "Types.h"

template <typename Field>
class MatrixSlice {
  static_assert(!std::is_reference_v<Field>);

  std::span<Field> data_;

  size_t rows_;
  size_t cols_;

  std::pair<size_t, size_t> rows_range_;
  std::pair<size_t, size_t> cols_range_;

  MatrixSlice(std::span<Field> data, size_t rows, size_t cols,
              std::pair<size_t, size_t> rows_range,
              std::pair<size_t, size_t> cols_range)
      : data_(data),
        rows_(rows),
        cols_(cols),
        rows_range_(rows_range),
        cols_range_(cols_range) {}

  static std::pair<size_t, size_t> as_range(
      const std::pair<size_t, size_t>& value) {
    return value;
  }

  static std::pair<size_t, size_t> as_range(size_t value) {
    return {value, value + 1};
  }

  size_t get_index(size_t row, size_t col) const {
    return col + cols_range_.first + cols_ * (row + rows_range_.first);
  }

  bool is_inside(size_t row, size_t col) const {
    return row + rows_range_.first < rows_range_.second &&
           col + cols_range_.first < cols_range_.second;
  }

  bool is_inside(std::pair<size_t, size_t> rows_range,
                 std::pair<size_t, size_t> cols_range) const {
    return rows_range.second + rows_range_.first <= rows_range_.second &&
           cols_range.second + cols_range_.first <= cols_range_.second;
  }

  template <typename U, typename V>
  decltype(auto) element_access_impl(U rows, V cols) {
    if constexpr (std::is_integral_v<U> && std::is_integral_v<V>) {
      if (!is_inside(rows, cols)) {
        throw std::invalid_argument("The index is outside of boundaries.");
      }

      return data_[get_index(rows, cols)];
    } else {
      std::pair<size_t, size_t> rows_range = as_range(rows);
      std::pair<size_t, size_t> cols_range = as_range(cols);

      if (!is_inside(rows_range, cols_range)) {
        throw std::invalid_argument("The index is outside of boundaries.");
      }

      rows_range.first += rows_range_.first;
      rows_range.second += rows_range_.first;

      cols_range.first += cols_range_.first;
      cols_range.second += cols_range_.first;

      return MatrixSlice(data_, rows_, cols_, rows_range, cols_range);
    }
  }

  MatrixSlice& apply_elementwise(MatrixSlice<const Field> other, auto&& functor)
    requires(!std::is_const_v<Field>)
  {
    if (other.shape() != shape()) {
      throw DimensionsException("Matrices must have same shape.");
    }

    auto [n, d] = shape();

    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < d; ++j) {
        functor((*this)[i, j], other[i, j]);
      }
    }

    return *this;
  }

 public:
  MatrixSlice& operator=(const MatrixSlice& other) {
    apply_elementwise(other,
                      [](Field& left, const Field& right) { left = right; });

    return *this;
  }

  MatrixSlice& operator=(const MatrixSlice<const Field>& other)
    requires(!std::is_const_v<Field>)
  {
    apply_elementwise(other,
                      [](Field& left, const Field& right) { left = right; });

    return *this;
  }

  Field& operator[](size_t row, size_t col) {
    return element_access_impl(row, col);
  }

  MatrixSlice operator[](std::pair<size_t, size_t> rows, size_t col) {
    return element_access_impl(rows, col);
  }

  MatrixSlice operator[](size_t row, std::pair<size_t, size_t> cols) {
    return element_access_impl(row, cols);
  }

  MatrixSlice operator[](std::pair<size_t, size_t> rows,
                         std::pair<size_t, size_t> cols) {
    return element_access_impl(rows, cols);
  }

  std::pair<size_t, size_t> shape() const {
    return {get_height(), get_width()};
  }

  size_t get_width() const { return cols_range_.second - cols_range_.first; }
  size_t get_height() const { return rows_range_.second - rows_range_.first; }

  // TODO: this can be replaced by carefully designed expression templates
  // https://en.wikipedia.org/wiki/Expression_templates
  MatrixSlice& add_mul(MatrixSlice<const Field> other, const Field& alpha) {
    return apply_elementwise(other, [&alpha](Field& left, const Field& right) {
      left += right * alpha;
    });
  }

  MatrixSlice& sub_mul(MatrixSlice<const Field> other, const Field& alpha) {
    return apply_elementwise(other, [&alpha](Field& left, const Field& right) {
      left -= right * alpha;
    });
  }

  MatrixSlice& operator+=(MatrixSlice<const Field> other) {
    return apply_elementwise(
        other, [](Field& left, const Field& right) { left += right; });
  }

  MatrixSlice& operator-=(MatrixSlice<const Field> other) {
    return apply_elementwise(
        other, [](Field& left, const Field& right) { left -= right; });
  }

  MatrixSlice& operator*=(const Field& value)
    requires(!std::is_const_v<Field>)
  {
    auto [n, d] = shape();

    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < d; ++j) {
        (*this)[i, j] *= value;
      }
    }

    return *this;
  }

  operator MatrixSlice<const Field>() const {
    return MatrixSlice<const Field>(data_, rows_, cols_, rows_range_,
                                    cols_range_);
  }

  friend Matrix<std::remove_const_t<Field>>;
  friend MatrixSlice<std::decay_t<Field>>;
};

namespace std {
template <typename Field>
void swap(MatrixSlice<Field> lhs, MatrixSlice<Field> rhs) {
  if (lhs.shape() != rhs.shape()) {
    throw DimensionsException("Matrices must have same shape for swap.");
  }

  auto [n, d] = lhs.shape();

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < d; ++j) {
      std::swap(lhs[i, j], rhs[i, j]);
    }
  }
}
}  // namespace std

template <typename Field>
class Matrix {
  std::vector<Field> data_;
  size_t rows_count_;
  size_t cols_count_;

  size_t get_index(size_t row, size_t col) const {
    return col + cols_count_ * row;
  }

 public:
  // Matrix(const Matrix& other)
  //     : data_(other.data_),
  //       rows_count_(other.rows_count_),
  //       cols_count_(other.cols_count_) {
  //   std::cout << "Matrix copy" << std::endl;
  // }
  //
  // Matrix(Matrix&& other)
  //     : data_(std::move(other.data_)),
  //       rows_count_(other.rows_count_),
  //       cols_count_(other.cols_count_) {
  //   std::cout << "Matrix move" << std::endl;
  // }
  //
  // Matrix& operator=(const Matrix& other) {
  //   std::cout << "Matrix copy assignment" << std::endl;
  //   data_ = other.data_;
  //   rows_count_ = other.rows_count_;
  //   cols_count_ = other.cols_count_;
  //
  //   return *this;
  // }
  //
  // Matrix& operator=(Matrix&& other) {
  //   std::cout << "Matrix move assignment" << std::endl;
  //   data_ = std::move(other.data_);
  //   rows_count_ = other.rows_count_;
  //   cols_count_ = other.cols_count_;
  //
  //   return *this;
  // }

  Matrix() : Matrix(0, 0) {}

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

  Matrix(MatrixSlice<Field> slice)
      : Matrix(slice.get_height(), slice.get_width()) {
    MatrixSlice<Field>(*this) = slice;
  }

  Matrix(MatrixSlice<const Field> slice)
      : Matrix(slice.get_height(), slice.get_width()) {
    MatrixSlice<Field>(*this) = slice;
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

  static Matrix item(Field value) { return Matrix(1, 1, std::move(value)); }

  template <typename Functor>
  static Matrix elementwise(MatrixSlice<const Field> left,
                            MatrixSlice<const Field> right, Functor functor) {
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

  operator MatrixSlice<Field>() {
    return MatrixSlice<Field>(data_, rows_count_, cols_count_, {0, rows_count_},
                              {0, cols_count_});
  }

  operator MatrixSlice<const Field>() const {
    return MatrixSlice<const Field>(data_, rows_count_, cols_count_,
                                    {0, rows_count_}, {0, cols_count_});
  }

  // TODO: maybe make this shorter
  Field& operator[](size_t row, size_t col) {
    return MatrixSlice<Field>(*this)[row, col];
  }

  MatrixSlice<Field> operator[](std::pair<size_t, size_t> rows, size_t col) {
    return MatrixSlice<Field>(*this)[rows, col];
  }

  MatrixSlice<Field> operator[](size_t row, std::pair<size_t, size_t> cols) {
    return MatrixSlice<Field>(*this)[row, cols];
  }

  MatrixSlice<Field> operator[](std::pair<size_t, size_t> rows,
                                std::pair<size_t, size_t> cols) {
    return MatrixSlice<Field>(*this)[rows, cols];
  }

  const Field& operator[](size_t row, size_t col) const {
    return MatrixSlice<const Field>(*this)[row, col];
  }

  MatrixSlice<const Field> operator[](std::pair<size_t, size_t> rows,
                                      size_t col) const {
    return MatrixSlice<const Field>(*this)[rows, col];
  }

  MatrixSlice<const Field> operator[](size_t row,
                                      std::pair<size_t, size_t> cols) const {
    return MatrixSlice<const Field>(*this)[row, cols];
  }

  MatrixSlice<const Field> operator[](std::pair<size_t, size_t> rows,
                                      std::pair<size_t, size_t> cols) const {
    return MatrixSlice<const Field>(*this)[rows, cols];
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
};

template <class F>
Matrix(MatrixSlice<F>) -> Matrix<std::remove_const_t<F>>;

template <MatrixLike L, MatrixLike R>
Matrix<common_field_t<L, R>> operator*(L&& left, R&& right) {
  auto [n1, m1] = left.shape();
  auto [n2, m2] = right.shape();

  if (m1 != n2) {
    throw std::invalid_argument("Arguments shapes are not compatible.");
  }

  auto result = Matrix<common_field_t<L, R>>(n1, m2);

  for (size_t i = 0; i < n1; ++i) {
    for (size_t j = 0; j < m2; ++j) {
      for (size_t k = 0; k < n2; ++k) {
        result[i, j] += left[i, k] * right[k, j];
      }
    }
  }

  return result;
}

template <MatrixLike L, MatrixLike R>
Matrix<common_field_t<L, R>> operator-(L&& left, R&& right) {
  using Field = common_field_t<L, R>;

  return Matrix<Field>::elementwise(
      left, right, [](const Field& l, const Field& r) { return l - r; });
}

template <MatrixLike L, MatrixLike R>
Matrix<common_field_t<L, R>> operator+(L&& left, R&& right) {
  using Field = common_field_t<L, R>;

  return Matrix<Field>::elementwise(
      left, right, [](const Field& l, const Field& r) { return l + r; });
}

std::ostream& operator<<(std::ostream& os, MatrixLike auto&& matrix) {
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
    os << "{";
    for (size_t j = 0; j < m; ++j) {
      os << std::right << std::setw(max_length) << results[i * m + j] << ", ";
    }

    os << "},\n";
  }

  return os;
}

namespace linalg {
template <MatrixLike Head, MatrixLike... Tail>
Matrix<common_field_t<Head, Tail...>> vstack(Head&& topmost, Tail&&... rest) {
  if (((topmost.get_width() != rest.get_width()) || ...)) {
    throw std::runtime_error(
        "Matrices must have equal width to stack them vertically.");
  }

  size_t width = topmost.get_width();
  Matrix<common_field_t<Head, Tail...>> result(
      (topmost.get_height() + ... + rest.get_height()), width);

  size_t row = 0;

  // not a snake!
  auto adder = [&result, &row, width](auto&& part) {
    for (size_t i = 0; i < part.get_height(); ++i) {
      for (size_t j = 0; j < width; ++j) {
        result[row + i, j] = part[i, j];
      }
    }

    row += part.get_height();
  };

  adder(topmost);
  (adder(rest), ...);

  return result;
}

template <MatrixLike Head, MatrixLike... Tail>
Matrix<common_field_t<Head, Tail...>> hstack(Head&& leftmost, Tail&&... rest) {
  if (((leftmost.get_height() != rest.get_height()) || ...)) {
    throw std::runtime_error(
        "Matrices must have equal height to stack them horizontally.");
  }

  size_t height = leftmost.get_height();
  Matrix<common_field_t<Head, Tail...>> result(
      height, (leftmost.get_width() + ... + rest.get_width()));

  size_t col = 0;

  // not a snake!
  auto adder = [&result, &col, height](auto&& part) {
    for (size_t i = 0; i < height; ++i) {
      for (size_t j = 0; j < part.get_width(); ++j) {
        result[i, col + j] = part[i, j];
      }
    }

    col += part.get_width();
  };

  adder(leftmost);
  (adder(rest), ...);

  return result;
}

template <MatrixLike T>
Matrix<matrix_field_t<T>> transposed(T&& matrix) {
  auto [n, d] = matrix.shape();

  Matrix<matrix_field_t<T>> result(d, n);

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < d; ++j) {
      result[j, i] = matrix[i, j];
    }
  }

  return result;
}

template <MatrixLike L, MatrixLike R>
common_field_t<L, R> dot(L&& left, R&& right) {
  size_t n = left.get_height();

  if (left.shape() != std::pair{n, 1} || right.shape() != std::pair{n, 1}) {
    throw DimensionsException(
        "Matrices must have shape (n, 1) for dot product.");
  }

  common_field_t<L, R> result;

  for (size_t i = 0; i < n; ++i) {
    result += left[i, 0] * right[i, 0];
  }

  return result;
}
}  // namespace linalg
