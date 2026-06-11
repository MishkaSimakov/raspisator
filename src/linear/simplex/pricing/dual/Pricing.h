#pragma once

#include <optional>

#include "linear/simplex/StateView.h"
#include "linear/simplex/Types.h"

namespace simplex {

template <typename Field>
class DualPricing {
 public:
  virtual std::optional<LeavingVariable> get_dual_leaving(
      StateView<Field> simplex) = 0;

  virtual ~DualPricing() = default;
};

}  // namespace simplex
