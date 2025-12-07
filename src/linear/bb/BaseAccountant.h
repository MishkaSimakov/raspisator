#pragma once

#include "linear/model/LP.h"
#include "linear/simplex/BaseAccountant.h"

template <typename Field>
struct BaseAccountant {
  void set_root(size_t id) {}

  void pruned_by_infeasibility(size_t id) {}
  void pruned_by_integrality(size_t id,
                             const FiniteLPSolution<Field>& solution) {}
  void pruned_by_bounds(size_t id, const FiniteLPSolution<Field>& solution) {}

  void set_branching_variable(size_t id,
                              const FiniteLPSolution<Field>& solution,
                              size_t variable_id) {}

  void set_left_child(size_t parent_id, size_t child_id) {}
  void set_right_child(size_t parent_id, size_t child_id) {}

  void simplex_run(size_t node_id, const SimplexResult<Field>&) {}
  void strong_branching_simplex_run(size_t node_id, size_t variable_id,
                                    const SimplexResult<Field>&) {}
};
