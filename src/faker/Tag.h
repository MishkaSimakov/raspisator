#pragma once

#include <cstdint>
#include <utility>

#include "faker/Instance.h"

namespace faker {

enum class Tag : uint64_t {
  // sizes
  TINY = 1 << 0,
  SMALL = 1 << 1,

  // solution type
  FEASIBLE = 1 << 2,
  INFEASIBLE = 1 << 3,
  UNBOUNDED = 1 << 4,

  // problem type
  LP = 1 << 5,
  MILP = 1 << 6,
  STANDARD_LP = 1 << 7,
  STANDARD_MILP = 1 << 8,

  // known properties
  KNOWN_OPTIMAL_OBJECTIVE = 1 << 9,

  // all variables have finite lower and upper bounds
  ALL_VARIABLES_BOUNDED = 1 << 10,
};

template <typename Field>
struct TaggedInstance {
  Tag tag;
  Instance<Field> instance;
};

constexpr Tag operator|(Tag a, Tag b) {
  return static_cast<Tag>(std::to_underlying(a) | std::to_underlying(b));
}

constexpr Tag operator&(Tag a, Tag b) {
  return static_cast<Tag>(std::to_underlying(a) & std::to_underlying(b));
}

constexpr bool has_tag(Tag mask, Tag flag) { return (mask & flag) != Tag{0}; }

}  // namespace faker
