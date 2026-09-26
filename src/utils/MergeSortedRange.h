#pragma once

#include <functional>
#include <span>
#include <tuple>

// Requires that x and y are sorted in strictly ascending order relative to
// Comparator.
template <typename T, typename Comparator = std::less<>>
class MergeSortedRange {
  using SparseVectorView = std::span<const std::pair<size_t, T>>;

  [[no_unique_address]]
  Comparator comparator_;

  SparseVectorView x_;
  SparseVectorView y_;

  class Iterator {
   public:
    using value_type = std::tuple<size_t, T, T>;
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::forward_iterator_tag;

    Iterator(SparseVectorView x, SparseVectorView y, size_t xi, size_t yi,
             Comparator comparator)
        : x_(x), y_(y), xi_(xi), yi_(yi), comparator_(comparator) {}

    value_type operator*() const {
      if (yi_ >= y_.size() ||
          (xi_ < x_.size() && comparator_(x_[xi_].first, y_[yi_].first))) {
        return {x_[xi_].first, x_[xi_].second, T()};
      }

      if (xi_ >= x_.size() || comparator_(y_[yi_].first, x_[xi_].first)) {
        return {y_[yi_].first, T(), y_[yi_].second};
      }

      return {x_[xi_].first, x_[xi_].second, y_[yi_].second};
    }

    Iterator& operator++() {
      if (yi_ >= y_.size() ||
          (xi_ < x_.size() && comparator_(x_[xi_].first, y_[yi_].first))) {
        ++xi_;
      } else if (xi_ >= x_.size() ||
                 comparator_(y_[yi_].first, x_[xi_].first)) {
        ++yi_;
      } else {
        ++xi_;
        ++yi_;
      }

      return *this;
    }

    bool operator==(const Iterator& other) const {
      return xi_ == other.xi_ && yi_ == other.yi_;
    }

    bool operator!=(const Iterator& other) const { return !(*this == other); }

   private:
    [[no_unique_address]]
    Comparator comparator_;

    SparseVectorView x_;
    SparseVectorView y_;

    size_t xi_;
    size_t yi_;
  };

 public:
  using iterator = Iterator;

  MergeSortedRange(SparseVectorView x, SparseVectorView y,
                   Comparator comparator = {})
      : x_(x), y_(y), comparator_(comparator) {}

  auto begin() const { return Iterator(x_, y_, 0, 0, comparator_); }

  auto end() const {
    return Iterator(x_, y_, x_.size(), y_.size(), comparator_);
  }
};
