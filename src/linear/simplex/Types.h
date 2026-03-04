#pragma once

#include "CyclingDetector.h"
#include "linear/matrix/Matrix.h"
#include "linear/model/LP.h"
#include "linear/sparse/LU.h"

namespace simplex {

template <typename Field>
struct IterationState {
  // Current iteration index
  size_t iteration_index;

  // Last iteration when cycling was detected, std::nullopt if cycling was not
  // detected yet.
  std::optional<size_t> last_cycling_iteration;

  // Current objective value
  Field objective;

  // Current basic variables
  std::vector<size_t> basic_variables;

  // Current variables states
  std::vector<VariableState> variables_states;

  // Current values of basic variables.
  // basic_point[i, 0] is the value of basic_variables[i]
  Matrix<Field> basic_point;

  // During simplex iterations it is guaranteed that this pointer is valid.
  const Bounds<Field>* bounds;

  // This LUPA instance stores the inverse of the matrix formed by taking basic
  // columns from A.
  linalg::LUPA<Field> lupa;

  CyclingDetector<Field> cycling;

  explicit IterationState(const CSCMatrix<Field>& A)
      : basic_variables(A.shape().first),
        variables_states(A.shape().second),
        basic_point(A.shape().first, 1),
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
