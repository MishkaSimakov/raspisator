#pragma once

#include "BaseView.h"
#include "OwningView.h"
#include "RefView.h"

#include "linalg/Concepts.h"

namespace linalg::detail {

template <MatrixRange M>
auto all(M&& matrix) {
  constexpr bool ref_suitable = requires(M matrix) { RefView(matrix); };

  if constexpr (std::is_base_of_v<BaseView, std::decay_t<M>>) {
    return std::forward<M>(matrix);
  } else if constexpr (ref_suitable) {
    return RefView(matrix);
  } else {
    return OwningView(matrix);
  }
}

template <MatrixRange M>
using all_t = decltype(all(std::declval<M>()));

}  // namespace linalg::detail
