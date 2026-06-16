#pragma once

#include <optional>
#include <string>

namespace simplex {

enum class Status {
  OPTIMAL,
  INFEASIBLE,
  UNBOUNDED,
  ITERATIONS_LIMIT,
};

template <typename Field>
struct Result {
  Status status;
  size_t iterations_count;
  std::optional<Field> objective;
};

inline std::string to_string(Status status) {
  switch (status) {
    case Status::OPTIMAL:
      return "OPTIMAL";
    case Status::INFEASIBLE:
      return "INFEASIBLE";
    case Status::UNBOUNDED:
      return "UNBOUNDED";
    case Status::ITERATIONS_LIMIT:
      return "ITERATIONS_LIMIT";
    default:
      std::unreachable();
  }
}

}  // namespace simplex
