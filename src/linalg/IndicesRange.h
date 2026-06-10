#pragma once

#include "Concepts.h"

namespace linalg {

// Returns view with values {from, from + 1, ..., to - 1}.
// Note: @to is excluded!
inline auto seq(size_t from, size_t to) {
  return std::ranges::iota_view(from, to);
}

}  // namespace linalg
