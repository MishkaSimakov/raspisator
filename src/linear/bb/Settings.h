#pragma once

#include <optional>

template <typename Field>
struct BranchAndBoundSettings {
  std::optional<size_t> max_nodes = std::nullopt;

  Field initial_pseudocost = 1;
  size_t reliability_parameter = 4;
  Field score_factor = static_cast<Field>(1) / 6;
};
