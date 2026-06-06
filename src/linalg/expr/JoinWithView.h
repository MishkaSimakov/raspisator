#pragma once

// This view is available in C++23 standard library, but my compiler doesn't
// support it, so I've implemented it on my own.
// This implementation may not fully cover the standard!
// Borrowed many things from llvm implementation.

#include <ranges>
#include <type_traits>
#include <variant>

namespace linalg::detail {

template <typename Left, typename Right>
class JoinWithView : std::ranges::view_interface<JoinWithView<Left, Right>> {
  Left left_;
  Right right_;

  template <bool IsConst>
  class Iterator;

  using LeftIter = std::ranges::iterator_t<Left>;
  using RightIter = std::ranges::iterator_t<Right>;

 public:
  JoinWithView(Left left, Right right)
      : left_(std::move(left)), right_(std::move(right)) {}

  Iterator<false> begin() {
    return Iterator<false>(std::variant<LeftIter, RightIter>(
        std::in_place_index<0>, std::ranges::begin(left_)));
  }

  Iterator<true> begin() const {
    return Iterator<true>(std::variant<LeftIter, RightIter>(
        std::in_place_index<0>, std::ranges::begin(left_)));
  }

  Iterator<false> end() {
    return Iterator<false>(std::variant<LeftIter, RightIter>(
        std::in_place_index<1>, std::ranges::end(right_)));
  }

  Iterator<true> end() const {
    return Iterator<true>(std::variant<LeftIter, RightIter>(
        std::in_place_index<1>, std::ranges::end(right_)));
  }
};

template <typename Left, typename Right>
template <bool IsConst>
class JoinWithView<Left, Right>::Iterator {
  using LeftIter = std::ranges::iterator_t<Left>;
  using RightIter = std::ranges::iterator_t<Right>;

  std::variant<LeftIter, RightIter> iter_;

  explicit Iterator(std::variant<LeftIter, RightIter> iter) : iter_(iter) {}

 public:
  using value_type = std::common_type_t<std::iter_value_t<LeftIter>,
                                        std::iter_value_t<RightIter>>;
  using difference_type = std::common_type_t<std::iter_difference_t<LeftIter>,
                                             std::iter_difference_t<RightIter>>;

  decltype(auto) operator*() const {
    return std::visit([](auto& it) { return *it; }, iter_);
  }

  Iterator& operator++() { return *this; }

  Iterator operator++(int) {
    auto copy = *this;
    ++*this;
    return copy;
  }

  bool operator==(const Iterator&) const = default;

  friend JoinWithView;
};

template <typename Left, typename Right>
JoinWithView(Left, Right)
    -> JoinWithView<std::views::all_t<Left>, std::views::all_t<Right>>;

}  // namespace linalg::detail
