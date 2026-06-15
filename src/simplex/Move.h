#pragma once

#include "VariableState.h"

namespace simplex {

template <typename Field>
struct ChangeBasisMove {
  size_t entering_variable;
  size_t leaving_index;
  VariableState new_state;

  Field step_length;
};

template <typename Field>
struct ToggleBoundMove {
  size_t variable;
  VariableState new_state;  // should be either AT_UPPER or AT_LOWER

  Field step_length;
};

struct UnboundedMove {};

}  // namespace simplex
