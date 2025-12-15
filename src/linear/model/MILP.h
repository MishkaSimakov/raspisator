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
using MILPSolution = std::variant<FiniteMILPSolution<Field>, NoFiniteSolution,
                                  ReachedNodesLimit>;
