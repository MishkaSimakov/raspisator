#pragma once

#include <variant>

#include "linear/matrix/Matrix.h"

enum class VariableType { INTEGER, SLACK, REAL };

// c x -> max, s.t. Ax = b, l <= x <= u
// and all integer_indices variables are integer
template <typename Field>
struct MILPProblem {
  Matrix<Field> A;
  Matrix<Field> b;
  Matrix<Field> c;

  std::vector<VariableType> variables_;
  std::vector<Field> lower_bounds;
  std::vector<Field> upper_bounds;
};

template <typename Field>
struct FiniteMILPSolution {
  Matrix<Field> point;
  Field value;
};

// TODO: for now MILP solver does not distinguish between NoFeasibleElements and
// InfiniteSolution
struct NoFiniteSolution {};

struct ReachedNodesLimit {};

template <typename Field>
using MILPSolution = std::variant<FiniteMILPSolution<Field>, NoFiniteSolution,
                                  ReachedNodesLimit>;
