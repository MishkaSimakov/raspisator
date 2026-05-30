#pragma once

#include <optional>

#include "linear/simplex/State.h"

namespace simplex {

// Pricing class is initialized each time simplex method starts, and it is kept
// until simplex has finished.
template <typename Field>
class PrimalPricing {
 public:
  virtual std::optional<size_t> get_primal_entering(
      detail::State<Field> simplex) = 0;

  virtual ~PrimalPricing() = default;
};

}  // namespace simplex
