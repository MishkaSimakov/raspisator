#pragma once

#include <format>
#include <iostream>

#include "Types.h"

namespace details {

struct NPYFilePrefix {
  unsigned char magic[6] = {0x93, 'N', 'U', 'M', 'P', 'Y'};
  unsigned char major = 0x01;
  unsigned char minor = 0x00;
};

template <typename Field>
struct NumpyFieldCode {};

template <>
struct NumpyFieldCode<double> : std::integral_constant<char, 'd'> {};

void write_header(std::ostream& os, MatrixLike auto&& matrix) {
  NPYFilePrefix prefix{};
  os.write(reinterpret_cast<char*>(&prefix), sizeof(prefix));

  auto field_code = NumpyFieldCode<matrix_field_t<decltype(matrix)>>::value;
  auto [n, d] = matrix.shape();

  std::string header = std::format(
      "{{'descr': '{}', 'fortran_order': False, 'shape': ({}, {})}}",
      field_code, n, d);

  uint16_t header_size = header.size();

  // pad header with whitespaces so that len(prefix) + 2 + len(header) % 64 == 0
  uint16_t total_size = header_size + 2 + sizeof(NPYFilePrefix);

  if (total_size % 64 != 0) {
    uint16_t padding = 64 - total_size % 64;
    header_size += padding;
    header += std::string(padding, ' ');
  }

  os.write(reinterpret_cast<char*>(&header_size), sizeof(header_size));
  os.write(header.c_str(), header_size);
}

void write_matrix(std::ostream& os, MatrixLike auto&& matrix)
  requires(std::same_as<matrix_field_t<decltype(matrix)>, double>)
{
  auto [n, d] = matrix.shape();

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < d; ++j) {
      os.write(reinterpret_cast<char*>(&matrix[i, j]),
               sizeof(matrix_field_t<decltype(matrix)>));
    }
  }
}

}  // namespace details

namespace linalg {

void to_npy(std::ostream& os, MatrixLike auto&& matrix) {
  details::write_header(os, matrix);

  details::write_matrix(os, matrix);
}

}  // namespace linalg
