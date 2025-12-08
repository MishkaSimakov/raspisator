#pragma once

#include <optional>

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
class Minimum {
  std::optional<Field> minimum_;

 public:
  Minimum() : minimum_(std::nullopt) {}

  void record(Field value) {
    if (!minimum_ ||
        FieldTraits<Field>::is_strictly_positive(*minimum_ - value)) {
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

template <typename Field>
class ArgMinimum {
  std::optional<std::pair<size_t, Field>> minimum_;

 public:
  ArgMinimum() : minimum_(std::nullopt) {}

  void record(size_t index, Field value) {
    if (!minimum_ ||
        FieldTraits<Field>::is_strictly_positive(minimum_->second - value)) {
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

template <typename Field>
class KahanSum {
  Field compensation_ = 0;
  Field sum_ = 0;

 public:
  void add(Field value) {
    double y = value - compensation_;
    double t = sum_ + y;

    compensation_ = (t - sum_) - y;
    sum_ = t;
  }

  Field sum() const { return sum_; }
};
