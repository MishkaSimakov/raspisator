#pragma once

#include <optional>

#include "linear/simplex/State.h"
#include "linear/simplex/Types.h"

namespace simplex {

template <typename Field>
class DualPricing {
 public:
  virtual std::optional<LeavingVariable> get_dual_leaving(
      detail::State<Field> simplex) = 0;

  virtual ~DualPricing() = default;
};

}  // namespace simplex
