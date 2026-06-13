#pragma once

#include <variant>
#include <vector>

#include "linalg/Linalg.h"
#include "simplex/VariableState.h"

namespace simplex {

template <typename Field>
struct FiniteLPSolution {
  Vector<Field> point;
  Field value;

  std::vector<VariableState> variables;

  std::vector<size_t> get_basic_variables() const {
    std::vector<size_t> result;

    for (size_t i = 0; i < variables.size(); ++i) {
      if (variables[i] == VariableState::BASIC) {
        result.push_back(i);
      }
    }

    return result;
  }
};

struct NoFeasibleElements {};

template <typename Field>
struct ReachedIterationsLimit {
  // dual objective value on the last iteration
  Field value;
};

struct Unbounded {};

template <typename Field>
struct SimplexResult {
  size_t iterations_count;

  std::variant<FiniteLPSolution<Field>, NoFeasibleElements,
               ReachedIterationsLimit<Field>, Unbounded>
      solution;

  bool is_feasible() const {
    return std::holds_alternative<FiniteLPSolution<Field>>(solution);
  }
};

}  // namespace simplex
