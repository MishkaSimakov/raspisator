#pragma once

#include <vector>

#include "linalg/Linalg.h"
#include "problem/Bound.h"
#include "simplex/Tolerance.h"
#include "simplex/VariableState.h"

namespace simplex {

// Simplex state DTO used for pricing and ratio tests
template <typename Field>
struct StateView {
  size_t iteration;

  Field objective;

  const Vector<Field>& basic_point;
  const std::vector<Bound<Field>>& bounds;
  const std::vector<VariableState>& states;
  const std::vector<size_t>& basic_vars;

  Tolerance<Field> tolerance;
};

}  // namespace simplex
