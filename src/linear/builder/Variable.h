#pragma once

#include "fwd.h"

template <typename Field>
class Variable {
  size_t id_;

  explicit Variable(size_t id) : id_(id) {}

  friend class ProblemBuilder<Field>;
  friend class Expression<Field>;
};
