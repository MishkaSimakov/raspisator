#pragma once

#include <optional>

#include "simplex/Move.h"
#include "simplex/StateView.h"

namespace simplex {

// Pricing class is initialized each time simplex method starts, and it is kept
// until simplex has finished.
template <typename Field>
class PrimalPricing {
 public:
  virtual void init(const problem::StandardLP<Field>& problem,
                    linalg::LUPA<Field>& lupa,
                    const std::vector<VariableState>& var_states,
                    const std::vector<size_t>& basic_vars) {}

  virtual std::optional<size_t> get_primal_entering(
      StateView<Field> simplex, const Vector<Field>& reduced_cost) = 0;

  virtual void post_refactorization(
      const problem::StandardLP<Field>& problem, linalg::LUPA<Field>& lupa,
      const std::vector<VariableState>& var_states,
      const std::vector<size_t>& basic_vars) {}

  virtual void move(ChangeBasisMove<Field> move, StateView<Field> simplex,
                    const Vector<Field>& pivot_row,
                    const Vector<Field>& pivot_col) {}
  virtual void move(ToggleBoundMove<Field> move, StateView<Field> simplex) {}

  virtual ~PrimalPricing() = default;
};

}  // namespace simplex
