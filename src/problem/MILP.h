#pragma once

#include "LP.h"

namespace problem {

template <typename Field>
struct MILP : LP<Field> {
  std::vector<bool> is_integer;
  std::vector<bool> implied_is_integer;

  // for debugging purposes, throws if problem is not correct
  void validate() const {
    LP<Field>::validate();

    const auto [n, d] = this->matrix.shape();

    if (is_integer.size() != d) {
      throw std::runtime_error("Wrong integrality vector size.");
    }

    if (implied_is_integer.size() != d) {
      throw std::runtime_error("Wrong implied integrality vector size.");
    }
  }
};

}  // namespace problem
