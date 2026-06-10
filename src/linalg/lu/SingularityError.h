#pragma once

#include <stdexcept>

namespace linalg {

struct SingularityError final : std::runtime_error {
  SingularityError()
      : std::runtime_error("Matrix is possibly singular in LU decomposition") {}
};

}