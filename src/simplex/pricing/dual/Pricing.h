#pragma once

#include <optional>

#include "simplex/StateView.h"

namespace simplex {

struct LeavingVariable {
  // Index of leaving variable in basic_variables array
  size_t index;
  VariableState new_state;
};

template <typename Field>
class DualPricing {
 public:
  virtual std::optional<LeavingVariable> get_dual_leaving(
      StateView<Field> simplex) = 0;

  virtual ~DualPricing() = default;
};

}  // namespace simplex
