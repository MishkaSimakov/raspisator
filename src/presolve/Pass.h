#pragma once

#include <cassert>
#include <vector>

#include "problem/MILP.h"

namespace presolve {

template <typename Field>
class Pass {
  bool was_applied_{false};

 protected:
  void register_apply() {
    if (was_applied_) {
      throw std::runtime_error("The same pass was applied twice.");
    }

    was_applied_ = true;
  }

 public:
  // must be called only once
  virtual problem::MILP<Field> apply(problem::MILP<Field> problem) = 0;

  // may be called many times
  virtual Vector<Field> inverse(Vector<Field> solution) const = 0;

  virtual ~Pass() = default;
};

}  // namespace presolve
