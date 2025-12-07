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

struct InfiniteSolution {};

struct NoFeasibleElements {};

template <typename Field>
using LPSolution =
    std::variant<FiniteLPSolution<Field>, InfiniteSolution, NoFeasibleElements>;

template <typename Field>
struct SimplexResult {
  size_t iterations_count;

  std::variant<FiniteLPSolution<Field>, NoFeasibleElements> solution;
};

// basic feasible solution
template <typename Field>
struct BFS {
  Matrix<Field> point;
  std::vector<size_t> basic_variables;

  static BFS construct_nondegenerate(Matrix<Field> point) {
    if (point.get_width() != 1) {
      throw std::invalid_argument("shape of the point must be (d, 1).");
    }

    std::vector<size_t> basic_variables;
    for (size_t i = 0; i < point.get_height(); ++i) {
      if (point[i, 0] != 0) {
        basic_variables.push_back(i);
      }
    }

    return BFS(point, basic_variables);
  }
};
