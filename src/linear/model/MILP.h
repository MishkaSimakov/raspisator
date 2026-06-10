#pragma once

#include <variant>

#include "linalg/Linalg.h"

template <typename Field>
struct FiniteMILPSolution {
  Vector<Field> point;
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
