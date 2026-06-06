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

  //
  size_t size() const { return this->rows(); }

  //
  void resize(size_t new_size) { this->resize(new_size, 1); }

  //
  using Matrix<Field>::operator[];

  Field& operator[](size_t row) { return (*this)[row, 0]; }
  const Field& operator[](size_t row) const { return (*this)[row, 0]; }
};

}  // namespace linalg
