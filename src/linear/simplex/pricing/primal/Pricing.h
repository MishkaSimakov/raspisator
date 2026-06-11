#pragma once

#include <optional>

#include "linear/simplex/StateView.h"

namespace simplex {

// Pricing class is initialized each time simplex method starts, and it is kept
// until simplex has finished.
template <typename Field>
class PrimalPricing {
 public:
  virtual std::optional<size_t> get_primal_entering(
      StateView<Field> simplex, const Vector<Field>& reduced_cost) = 0;

  virtual ~PrimalPricing() = default;
};

}  // namespace simplex
