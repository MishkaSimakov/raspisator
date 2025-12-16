#pragma once
#include "linear/problem/MILPProblem.h"

template <typename Field>
class BaseOptimizer {
 public:
  virtual MILPProblem<Field> apply(MILPProblem<Field> problem) = 0;

  virtual Matrix<Field> inverse(const Matrix<Field>& point) = 0;

  virtual ~BaseOptimizer() = default;
};
