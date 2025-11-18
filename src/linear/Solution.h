#pragma once

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

template<typename Field>
struct FiniteMILPSolution {
  Matrix<Field> point;
  Field value;
};