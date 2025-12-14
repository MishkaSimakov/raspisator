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

  std::vector<VariableType> variables;
  std::vector<Field> lower_bounds;
  std::vector<Field> upper_bounds;

  MILPProblem(const Matrix<Field>& A, const Matrix<Field>& b,
              const Matrix<Field>& c,
              const std::vector<VariableType>& variables,
              const std::vector<Field>& lower_bounds,
              const std::vector<Field>& upper_bounds)
      : A(A),
        b(b),
        c(c),
        variables(variables),
        lower_bounds(lower_bounds),
        upper_bounds(upper_bounds) {
    // check dimensions
    auto [n, d] = A.shape();

    if (b.shape() != std::pair{n, 1}) {
      throw DimensionsException("Wrong b shape.");
    }

    if (c.shape() != std::pair{1, d}) {
      throw DimensionsException("Wrong c shape.");
    }

    if (variables.size() != d) {
      throw DimensionsException("Wrong variables size.");
    }

    if (lower_bounds.size() != d) {
      throw DimensionsException("Wrong lower_bounds size.");
    }

    if (upper_bounds.size() != d) {
      throw DimensionsException("Wrong upper_bounds size.");
    }
  }
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
