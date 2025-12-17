#pragma once

#include <variant>

#include "linear/matrix/Matrix.h"

template <typename Field>
struct FiniteMILPSolution {
  Matrix<Field> point;
  Field value;
};

struct NoFiniteSolution {};

struct ReachedNodesLimit {};

template <typename Field>
struct BBRunResult {
  std::variant<FiniteMILPSolution<Field>, NoFiniteSolution, ReachedNodesLimit>
      solution;

  size_t nodes_count;
  double average_simplex_iterations;
};
