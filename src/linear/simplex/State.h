#pragma once

#include <vector>

#include "linalg/Linalg.h"
#include "linear/model/Bound.h"
#include "linear/model/LP.h"
#include "linear/simplex/Tolerance.h"

namespace simplex::detail {

// Simplex state DTO used for pricing and ratio tests
template <typename Field>
struct State {
  size_t iteration;

  Field objective;

  const Vector<Field>& basic_point;
  const std::vector<Bound<Field>>& bounds;
  const Vector<Field>& reduced_cost;
  const std::vector<VariableState>& states;
  const std::vector<size_t>& basic_vars;

  Tolerance<Field> tolerance;
};

}  // namespace simplex::detail
