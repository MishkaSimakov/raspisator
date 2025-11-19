#pragma once

#include <variant>

#include "linear/Matrix.h"

template <typename Field>
struct FiniteMILPSolution {
  Matrix<Field> point;
  Field value;
};

// TODO: for now MILP solver does not distinguish between NoFeasibleElements and
// InfiniteSolution
struct NoFiniteSolution {};

template <typename Field>
using MILPSolution = std::variant<FiniteMILPSolution<Field>, NoFiniteSolution>;
