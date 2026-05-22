#pragma once

#include <string>
#include <vector>

#include "linear/model/Bound.h"
#include "linear/sparse/CSCMatrix.h"

namespace problem {

template <typename Field>
struct LP {
  CSCMatrix<Field> matrix;
  std::vector<Field> cost;

  Field cost_offset;

  std::vector<Bound<Field>> var_bounds;
  std::vector<Bound<Field>> rhs_bounds;

  std::vector<Bound<Field>> implied_var_bounds;

  std::string name;
  std::string cost_name;

  std::vector<std::string> var_names;
  std::vector<std::string> row_names;
};

}  // namespace problem
