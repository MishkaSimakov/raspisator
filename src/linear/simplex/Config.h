#pragma once

#include <memory>
#include <optional>

#include "Accountant.h"
#include "Tolerance.h"
#include "linear/simplex/pricing/dual/Pricing.h"
#include "linear/simplex/pricing/primal/Pricing.h"
#include "pricing/dual/Dantzig.h"
#include "pricing/primal/MostInfeasible.h"

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

  std::unique_ptr<Accountant<Field>> accountant;

  Config()
      : primal_pricing(std::make_unique<PrimalMostInfeasible<Field>>()),
        dual_pricing(std::make_unique<DualDantzigPricing<Field>>()) {}

  Config& set_max_iterations(std::optional<size_t> max_iterations) & {
    this->max_iterations = max_iterations;
    return *this;
  }
  Config&& set_max_iterations(std::optional<size_t> max_iterations) && {
    this->max_iterations = max_iterations;
    return std::move(*this);
  }

  Config& set_validate_input(bool validate_input) & {
    this->validate_input = validate_input;
    return *this;
  }
  Config&& set_validate_input(bool validate_input) && {
    this->validate_input = validate_input;
    return std::move(*this);
  }

  template <typename T, typename... Args>
  Config& set_accountant(Args&&... args) & {
    this->accountant = std::make_unique<T>(std::forward<Args>(args)...);
    return *this;
  }
  template <typename T, typename... Args>
  Config&& set_accountant(Args&&... args) && {
    this->accountant = std::make_unique<T>(std::forward<Args>(args)...);
    return std::move(*this);
  }
};

}  // namespace simplex
