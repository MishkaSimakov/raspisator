#pragma once

#include <variant>
#include <vector>

#include "Bound.h"
#include "linear/matrix/Matrix.h"

enum class VariableState { AT_LOWER, AT_UPPER, BASIC };

template <typename Field>
class Bounds {
  std::vector<Bound<Field>> variables_bounds_;

 public:
  explicit Bounds(size_t size) : variables_bounds_(size) {}

  Bounds(const std::vector<Field>& lower, const std::vector<Field>& upper)
      : Bounds(lower.size()) {
    assert(lower.size() == upper.size());

    for (size_t i = 0; i < lower.size(); ++i) {
      variables_bounds_[i].lower = lower[i];
      variables_bounds_[i].upper = upper[i];
    }
  }

  Bound<Field>& operator[](size_t index) { return variables_bounds_[index]; }

  const Bound<Field>& operator[](size_t index) const {
    return variables_bounds_[index];
  }
};

template <typename Field>
struct FiniteLPSolution {
  Matrix<Field> point;
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
