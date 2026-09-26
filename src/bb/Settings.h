#pragma once

#include <optional>

template <typename Field>
struct BranchAndBoundSettings {
  std::optional<size_t> max_nodes = std::nullopt;

  std::optional<size_t> strong_branching_max_iterations_factor = 5;
  size_t strong_branching_min_iterations_limit = 0;

  // for reliability branching
  Field initial_pseudocost = 1;
  size_t reliability_parameter = 4;

  Field score_factor = static_cast<Field>(1) / 6;
};
