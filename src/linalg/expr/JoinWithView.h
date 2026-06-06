#pragma once

// This view is available in C++23 standard library, but my compiler doesn't
// support it, so I've implemented it on my own.
// This implementation may not fully cover the standard!
// Borrowed many things from llvm implementation.

#include <ranges>
#include <type_traits>
#include <variant>

namespace linalg::detail {

template <typename L, typename R>
class JoinWithView : std::ranges::view_interface<JoinWithView<L, R>> {
  L left_;
  R right_;

  template <bool IsConst>
  class Iterator;

  template <bool IsConst>
  class Sentinel;

  using LeftIter = std::ranges::iterator_t<L>;
  using RightIter = std::ranges::iterator_t<R>;

 public:
  JoinWithView(L left, R right)
      : left_(std::move(left)), right_(std::move(right)) {}

  Iterator<false> begin() {
    return Iterator<false>(
        std::variant<LeftIter, RightIter>(std::in_place_index<0>,
                                          std::ranges::begin(left_)),
        this);
  }

  Iterator<true> begin() const {
    return Iterator<true>(
        std::variant<LeftIter, RightIter>(std::in_place_index<0>,
                                          std::ranges::begin(left_)),
        this);
  }

  auto end() {
    if constexpr (std::ranges::common_range<R>) {
      return Iterator<false>(
          std::variant<LeftIter, RightIter>(std::in_place_index<1>,
                                            std::ranges::end(right_)),
          this);
    } else {
      return Sentinel<false>{std::ranges::end(right_)};
    }
  }

  auto end() const {
    if constexpr (std::ranges::common_range<R>) {
      return Iterator<true>(
          std::variant<LeftIter, RightIter>(std::in_place_index<1>,
                                            std::ranges::end(right_)),
          this);
    } else {
      return Sentinel<true>{std::ranges::end(right_)};
    }
  }
};

template <typename Left, typename Right>
template <bool IsConst>
class JoinWithView<Left, Right>::Iterator {
  using LeftIter = std::ranges::iterator_t<Left>;
  using RightIter = std::ranges::iterator_t<Right>;

  using Parent = std::conditional_t<IsConst, const JoinWithView, JoinWithView>;

  std::variant<LeftIter, RightIter> iter_;
  Parent* parent_;

  explicit Iterator(std::variant<LeftIter, RightIter> iter, Parent* parent)
      : iter_(iter), parent_(parent) {
    normalize();
  }

  void normalize() {
    if (iter_.index() == 0 &&
        std::get<0>(iter_) == std::ranges::end(parent_->left_)) {
      iter_.template emplace<1>(std::ranges::begin(parent_->right_));
    }
  }

 public:
  using value_type = std::common_type_t<std::iter_value_t<LeftIter>,
                                        std::iter_value_t<RightIter>>;
  using difference_type = std::common_type_t<std::iter_difference_t<LeftIter>,
                                             std::iter_difference_t<RightIter>>;

  Iterator() = default;

  decltype(auto) operator*() const {
    return std::visit([](auto& it) { return *it; }, iter_);
  }

  Iterator& operator++() {
    if (iter_.index() == 0) {
      ++std::get<0>(iter_);
      normalize();
    } else {
      ++std::get<1>(iter_);
    }

    return *this;
  }

  Iterator operator++(int) {
    auto copy = *this;
    ++*this;
    return copy;
  }

  bool operator==(const Iterator&) const = default;

  template <bool OtherIsConst>
  friend class Sentinel;

  friend JoinWithView;
};

template <typename Left, typename Right>
template <bool IsConst>
class JoinWithView<Left, Right>::Sentinel {
  using RightSentinel = std::ranges::sentinel_t<Right>;

  RightSentinel sentinel_;

  explicit Sentinel(RightSentinel sentinel) : sentinel_(std::move(sentinel)) {}

 public:
  Sentinel() = default;

  template <bool OtherIsConst>
  bool operator==(const Iterator<OtherIsConst>& it) const {
    return it.iter_.index() == 1 && std::get<1>(it.iter_) == sentinel_;
  }

  friend JoinWithView;
};

template <typename L, typename R>
JoinWithView(L&&, R&&)
    -> JoinWithView<std::views::all_t<L>, std::views::all_t<R>>;

}  // namespace linalg::detail
