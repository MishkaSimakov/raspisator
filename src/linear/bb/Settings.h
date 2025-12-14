#pragma once

#include <optional>

enum class PerturbationMode { DISABLED, CONSTANT, FOR_INTEGER_SOLUTION };

template <typename Field>
struct BranchAndBoundSettings {
  std::optional<size_t> max_nodes = std::nullopt;

  Field initial_pseudocost = 1;
  size_t reliability_parameter = 4;
  Field score_factor = static_cast<Field>(1) / 6;

  PerturbationMode perturbation = PerturbationMode::DISABLED;
  Field perturbation_value = static_cast<Field>(1e-5);
};
