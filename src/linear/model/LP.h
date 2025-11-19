#pragma once

#include <concepts>
#include <variant>
#include <vector>

#include "linear/Matrix.h"

template <typename Field>
struct FiniteLPSolution {
  Matrix<Field> point;
  std::vector<size_t> basic_variables;
  Field value;

  Matrix<Field> tableau;
};

struct InfiniteSolution {};

struct NoFeasibleElements {};

template <typename Field>
using LPSolution =
    std::variant<FiniteLPSolution<Field>, InfiniteSolution, NoFeasibleElements>;

template <typename T, typename Field>
concept LPSolver = requires(T solver) {
  requires std::constructible_from<T, Matrix<Field>, Matrix<Field>,
                                   Matrix<Field>>;
  { solver.solve() } -> std::same_as<LPSolution<Field>>;
};
