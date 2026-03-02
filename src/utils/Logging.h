#pragma once

#include <fstream>
#include <string_view>

#include "Paths.h"
#include "linear/matrix/Matrix.h"
#include "linear/matrix/NPY.h"
#include "linear/sparse/CSCMatrix.h"

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
void log_density(const CSCMatrix<Field>& matrix, std::string_view filename) {
  get_log_fstream(filename)
      << matrix.density() << ", " << matrix.nonzero_count() << "\n";
}

template <typename Field>
void log_density(const Matrix<Field>& matrix, std::string_view filename) {
  get_log_fstream(filename)
      << matrix.density() << ", " << matrix.nonzero_count() << "\n";
}

template <typename T>
void log_value(const T& value, std::string_view filename) {
  get_log_fstream(filename) << value << "\n";
}

template <typename Field>
void log_npy(const Matrix<Field>& matrix, std::string_view filename) {
  auto os = get_log_fstream(filename);

  linalg::to_npy(os, matrix);
}

}  // namespace logging
