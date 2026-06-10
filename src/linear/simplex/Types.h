#pragma once

#include "CyclingDetector.h"
#include "linalg/Linalg.h"
#include "linalg/lu/LUPA.h"
#include "linear/model/LP.h"

namespace simplex {

template <typename Field>
struct IterationState {
  // Current iteration index
  size_t iteration_index;

  // Current objective value
  Field objective;

  // Current basic variables
  std::vector<size_t> basic_variables;

  // Current variables states
  std::vector<VariableState> variables_states;

  // Current values of basic variables.
  // basic_point[i, 0] is the value of basic_variables[i]
  Vector<Field> basic_point;

  Vector<Field> reduced_cost;

  // During simplex iterations it is guaranteed that this pointer is valid.
  const Bounds<Field>* bounds;

  // This LUPA instance stores the inverse of the matrix formed by taking basic
  // columns from A.
  linalg::LUPA<Field> lupa;

  explicit IterationState(const linalg::CSCMatrix<Field>& A)
      : basic_variables(A.rows()),
        variables_states(A.cols()),
        basic_point(A.rows()),
        bounds(nullptr),
        lupa(A) {}

  std::pair<size_t, size_t> problem_shape() const {
    return {basic_variables.size(), variables_states.size()};
  }
};

struct LeavingVariable {
  // Index of leaving variable in basic_variables array
  size_t index;
  VariableState new_state;
};

}  // namespace simplex
