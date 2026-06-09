#pragma once

#include "Matrix.h"

#include "Arithmetics.h"

namespace linalg {

template <typename Field>
class Vector : public Matrix<Field> {
 public:
  using FieldType = Field;

  Vector() = default;

  explicit Vector(size_t rows) : Matrix<Field>(rows, 1) {}

  Vector(std::initializer_list<Field> il)
      : Matrix<Field>(Matrix<Field>::uninitialized(il.size(), 1)) {
    size_t row = 0;

    for (Field value : il) {
      (*this)[row, 0] = value;
      ++row;
    }
  }

  template <MatrixRange T>
    requires std::same_as<MatrixFieldType<T>, Field>
  Vector(T&& other) : Matrix<Field>(std::forward<T>(other)) {}

  //
  size_t size() const { return this->rows(); }

  //
  void resize(size_t new_size) { Matrix<Field>::resize(new_size, 1); }

  //
  using Matrix<Field>::operator[];

  Field& operator[](size_t row) { return (*this)[row, 0]; }
  const Field& operator[](size_t row) const { return (*this)[row, 0]; }

  template <IndicesRange RowRange>
  auto operator[](RowRange&& rows) const {
    return (*this)[std::forward<RowRange>(rows), std::views::single(size_t{0})];
  }

  template <typename F>
  void row_entries(size_t row, F&& f) const {
    f(0, (*this)[row, 0]);
  }

  template <typename F>
  void entries(F&& f) const {
    for (size_t row = 0; row < this->rows(); ++row) {
      f(row, 0, (*this)[row, 0]);
    }
  }
};

template <MatrixRange T>
Vector(T&&) -> Vector<MatrixFieldType<T>>;

}  // namespace linalg
