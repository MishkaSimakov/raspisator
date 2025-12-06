#pragma once

#include <optional>

template <typename Field>
struct BranchAndBoundSettings {
  std::optional<size_t> max_nodes = std::nullopt;
};
