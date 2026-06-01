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

  // When validate_input is true, then simplex will perform additional checks on
  // user input. The checks may be not exhaustive.
  // Input validation leads to worse performance.
  bool validate_input{false};

  Tolerance<Field> tolerance = kDefaultTolerance<Field>;

  std::unique_ptr<PrimalPricing<Field>> primal_pricing;
  std::unique_ptr<DualPricing<Field>> dual_pricing;
};

}  // namespace simplex
