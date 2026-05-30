#pragma once

#include <memory>
#include <optional>

#include "Tolerance.h"
#include "linear/simplex/pricing/dual/Pricing.h"
#include "linear/simplex/pricing/primal/Pricing.h"

namespace simplex {

template <typename Field>
struct Config {
  std::optional<size_t> max_iterations;

  // In strict mode simplex performs various checks of correctness. They lead to
  // worse performance.
  bool is_strict{false};

  Tolerance<Field> tolerance = kDefaultTolerance<Field>;

  std::unique_ptr<PrimalPricing<Field>> primal_pricing;
  std::unique_ptr<DualPricing<Field>> dual_pricing;
};

}  // namespace simplex
