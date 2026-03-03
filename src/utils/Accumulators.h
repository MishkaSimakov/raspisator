#pragma once

#include <optional>

#include "linear/FieldTraits.h"

template <typename Field>
class ArithmeticMean {
  Field sum_;
  size_t count_;

 public:
  ArithmeticMean() : sum_(0), count_(0) {}

  size_t count() const { return count_; }

  Field sum() const { return sum_; }

  Field mean() const { return sum_ / static_cast<Field>(count_); }

  void record(Field value) {
    sum_ += value;
    ++count_;
  }
};

template <typename Field>
class GeometricAverage {
  Field product_;
  size_t count_;

 public:
  GeometricAverage() : product_(1), count_(0) {}

  size_t count() const { return count_; }

  Field product() const { return product_; }

  Field average() const {
    if (count_ == 0) {
      return 1;
    }

    return std::pow(product_, 1. / count_);
  }

  void record(Field value) {
    product_ *= value;
    ++count_;
  }
};

template <typename Field>
struct FieldTraitsComparator {
  bool operator()(Field left, Field right) const {
    return FieldTraits<Field>::is_strictly_negative(left - right);
  }
};

template <typename Field, typename Comparator = FieldTraitsComparator<Field>>
class Minimum {
  std::optional<Field> minimum_;

  [[no_unique_address]]
  Comparator comparator_;

 public:
  Minimum() : minimum_(std::nullopt) {}

  void reset() { minimum_ = std::nullopt; }

  void record(Field value) {
    if (!minimum_ || comparator_(value, *minimum_)) {
      minimum_ = value;
    }
  }

  void record(std::optional<Field> value) {
    if (value) {
      record(*value);
    }
  }

  std::optional<Field> min() const { return minimum_; }
};

template <typename Field, typename Comparator = FieldTraitsComparator<Field>>
class Maximum {
  std::optional<Field> maximum_;

  [[no_unique_address]]
  Comparator comparator_;

 public:
  Maximum() : maximum_(std::nullopt) {}

  void reset() { maximum_ = std::nullopt; }

  void record(Field value) {
    if (!maximum_ || comparator_(*maximum_, value)) {
      maximum_ = value;
    }
  }

  void record(std::optional<Field> value) {
    if (value) {
      record(*value);
    }
  }

  std::optional<Field> max() const { return maximum_; }
};

// Calculates minimum i, s.t. a_i = min_j a_j
template <typename Field, typename Comparator = FieldTraitsComparator<Field>>
class ArgMinimum {
  std::optional<std::pair<size_t, Field>> minimum_;

  [[no_unique_address]]
  Comparator comparator_;

 public:
  ArgMinimum() : minimum_(std::nullopt) {}

  void record(size_t index, Field value) {
    if (!minimum_) {
      minimum_ = {index, value};
      return;
    }

    if (comparator_(value, minimum_->second)) {
      minimum_ = {index, value};
    } else if (!comparator_(minimum_->second, value) &&
               index < minimum_->first) {
      minimum_ = {index, value};
    }
  }

  void record(size_t index, std::optional<Field> value) {
    if (value) {
      record(index, *value);
    }
  }

  std::optional<Field> min() const {
    return minimum_.transform(
        [](std::pair<size_t, Field> value) { return value.second; });
  }

  std::optional<size_t> argmin() const {
    return minimum_.transform(
        [](std::pair<size_t, Field> value) { return value.first; });
  }
};

// Calculates minimum i, s.t. a_i = max_j a_j
template <typename Field, typename Comparator = FieldTraitsComparator<Field>>
class ArgMaximum {
  std::optional<std::pair<size_t, Field>> maximum_;

  [[no_unique_address]]
  Comparator comparator_;

 public:
  ArgMaximum() : maximum_(std::nullopt) {}

  void record(size_t index, Field value) {
    if (!maximum_) {
      maximum_ = {index, value};
      return;
    }

    if (comparator_(maximum_->second, value)) {
      maximum_ = {index, value};
    } else if (!comparator_(value, maximum_->second) &&
               index < maximum_->first) {
      maximum_ = {index, value};
    }
  }

  void record(size_t index, std::optional<Field> value) {
    if (value) {
      record(index, *value);
    }
  }

  std::optional<Field> max() const {
    return maximum_.transform(
        [](std::pair<size_t, Field> value) { return value.second; });
  }

  std::optional<size_t> argmax() const {
    return maximum_.transform(
        [](std::pair<size_t, Field> value) { return value.first; });
  }
};

template <typename Field>
class KahanSum {
  Field compensation_ = 0;
  Field sum_ = 0;

 public:
  void add(Field value) {
    Field y = value - compensation_;
    Field t = sum_ + y;

    compensation_ = (t - sum_) - y;
    sum_ = t;
  }

  Field sum() const { return sum_; }
};
