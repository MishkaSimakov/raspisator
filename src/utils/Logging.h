#pragma once

#include <fstream>

#include "linear/sparse/CSCMatrix.h"
#include "linear/matrix/Matrix.h"

namespace logging {

template<typename Field>
void log_density(const CSCMatrix<Field>& matrix, std::string_view filename) {
  std::ofstream os(filename, std::ofstream::app);

  os << matrix.density() << ", " << matrix.nonzero_count() << "\n";
}

template<typename Field>
void log_density(const Matrix<Field>& matrix, std::string_view filename) {
  std::ofstream os(filename, std::ofstream::app);

  os << matrix.density() << ", " << matrix.nonzero_count() << "\n";
}

}
