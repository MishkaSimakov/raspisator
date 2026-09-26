#pragma once

#include <vector>

#include "linalg/Linalg.h"
#include "linalg/lu/LUPA.h"
#include "problem/Bound.h"
#include "problem/StandardLP.h"
#include "simplex/Tolerance.h"
#include "simplex/VariableState.h"

namespace simplex {

// Simplex state DTO used for pricing and ratio tests
template <typename Field>
struct StateView {
  const problem::StandardLP<Field>& problem;

  size_t iteration;

  // Sometimes simplex intentionally repeat iteration for the same basis, and
  // this can confuse cycling detectors.
  bool intentional_repeat;

  Field objective;

  linalg::LUPA<Field>& lupa;

  const Vector<Field>& basic_point;
  const std::vector<VariableState>& states;
  const std::vector<size_t>& basic_vars;

  const Vector<Field>& reduced_cost;

  Tolerance<Field> tolerance;
};

}  // namespace simplex
