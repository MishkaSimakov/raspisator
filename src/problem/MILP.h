#pragma once

#include "LP.h"

namespace problem {

template <typename Field>
struct MILP : LP<Field> {
  std::vector<bool> is_integer;
  std::vector<bool> implied_is_integer;
};

}  // namespace problem
