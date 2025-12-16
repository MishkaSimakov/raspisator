#pragma once

#include <concepts>
#include <optional>
#include <variant>
#include <vector>

#include "linear/matrix/Matrix.h"
#include "linear/sparse/CSCMatrix.h"

enum class VariableState { AT_LOWER, AT_UPPER, BASIC };

template <typename Field>
struct FiniteLPSolution {
  Matrix<Field> point;
  std::vector<size_t> basic_variables;
  Field value;

  std::vector<VariableState> variables;
};

struct NoFeasibleElements {};

template <typename Field>
struct ReachedIterationsLimit {
  // dual objective value on the last iteration
  Field value;
};

template <typename Field>
struct SimplexResult {
  size_t iterations_count;

  std::variant<FiniteLPSolution<Field>, NoFeasibleElements,
               ReachedIterationsLimit<Field>>
      solution;

  bool is_feasible() const {
    return std::holds_alternative<FiniteLPSolution<Field>>(solution);
  }
};
