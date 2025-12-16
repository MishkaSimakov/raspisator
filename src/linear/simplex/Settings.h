#pragma once
#include <optional>

namespace simplex {

template <typename Field>
struct Settings {
  std::optional<size_t> max_iterations;
};

}  // namespace simplex
