#pragma once

#include <optional>

#include "linalg/Linalg.h"
#include "problem/MILP.h"

namespace faker {

enum class SolutionType {
  FEASIBLE,
  UNBOUNDED,
  INFEASIBLE,
  UNKNOWN,
};

enum class ProblemType {
  LP,
  StandardLP,
  MILP,
  StandardMILP,
};

template <typename Field>
struct Instance {
  problem::MILP<Field> problem;
  SolutionType solution_type;
  ProblemType problem_type;

  bool has_linearly_dependent_rows;

  // std::nullopt if unknown or absent
  // Note: may be unknown even if solution_type is FEASIBLE
  std::optional<Vector<Field>> feasible_point;
  std::optional<Field> optimal_objective;
};

}  // namespace faker
