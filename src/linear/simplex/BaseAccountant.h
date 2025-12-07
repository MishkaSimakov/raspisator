#pragma once

namespace simplex {

template <typename Field>
struct BaseAccountant {
  void store_iteration() {}

  void new_problem() {}
};

}  // namespace simplex
