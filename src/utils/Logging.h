#pragma once

#include <fstream>
#include <string_view>

#include "Paths.h"
#include "linalg/Linalg.h"
#include "linalg/NPY.h"

namespace logging {

inline std::ofstream get_log_fstream(std::string_view filename) {
  auto path = paths::log(filename);

  std::filesystem::create_directories(path.parent_path());

  std::ofstream os(path, std::ofstream::app);

  if (!os) {
    throw std::runtime_error("Failed to write to log file.");
  }

  return os;
}

template <typename Field>
void density(const CSCMatrix<Field>& matrix, std::string_view filename) {
  get_log_fstream(filename)
      << matrix.density() << ", " << matrix.nonzero_count() << "\n";
}

template <typename Field>
void density(const Matrix<Field>& matrix, std::string_view filename) {
  get_log_fstream(filename)
      << matrix.density() << ", " << matrix.nonzero_count() << "\n";
}

template <typename T>
void value(const T& value, std::string_view filename) {
  get_log_fstream(filename) << value << "\n";
}

template <typename Field>
void npy(const Matrix<Field>& matrix, std::string_view filename) {
  auto os = get_log_fstream(filename);

  linalg::to_npy(os, matrix);
}

inline void string(std::string_view text, std::string_view filename) {
  auto os = get_log_fstream(filename);
  os << text;
}

}  // namespace logging
