#pragma once

#include "BaseView.h"
#include "OwningView.h"
#include "RefView.h"

#include "linalg/Concepts.h"

namespace linalg::detail {

template <MatrixRange M>
auto all(M&& matrix) {
  if constexpr (std::is_base_of_v<BaseView, std::decay_t<M>>) {
    return std::forward<M>(matrix);
  } else if constexpr (std::is_reference_v<M>) {
    return RefView(std::forward<M>(matrix));
  } else {
    return OwningView(std::forward<M>(matrix));
  }
}

template <MatrixRange M>
using all_t = decltype(all(std::declval<M>()));

}  // namespace linalg::detail
