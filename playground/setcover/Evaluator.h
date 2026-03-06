#pragma once

#include "Types.h"

namespace setcover {

EvaluationResult evaluate(const Problem& problem, const Solution& solution) {
  // check that all elements are covered
  for (size_t i = 0; i < problem.elements_count; ++i) {
    bool is_covered = false;

    for (size_t set_index : solution.chosen_sets) {
      if (problem.sets[set_index].elements.contains(i)) {
        is_covered = true;
        break;
      }
    }

    if (!is_covered) {
      return EvaluationResult{
          .score = 0,
          .is_valid = false,
      };
    }
  }

  // calculate score
  size_t score = 0;
  for (size_t set_index : solution.chosen_sets) {
    score += problem.sets[set_index].cost;
  }

  return EvaluationResult{
      .score = score,
      .is_valid = true,
  };
}

}  // namespace setcover
