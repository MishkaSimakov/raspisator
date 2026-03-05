#pragma once
#include <optional>

namespace simplex {

template <typename Field>
struct Settings {
  std::optional<size_t> max_iterations;

  // In strict mode simplex performs various checks of correctness. They lead to
  // worse performance.
  bool is_strict{false};
};

}  // namespace simplex
