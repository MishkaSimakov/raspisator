#pragma once

#include <cassert>
#include <cmath>
#include <optional>

template <typename Field>
class ArithmeticMean {
  Field sum_;
  size_t count_;

 public:
  ArithmeticMean() : sum_(0), count_(0) {}

  void record(Field value) {
    sum_ += value;
    ++count_;
  }

  std::optional<Field> get() const {
    if (count_ == 0) {
      return std::nullopt;
    }

    return sum_ / count_;
  }
  bool has_value() const { return count_ != 0; }

  Field operator*() const {
    assert(count_ > 0);
    return sum_ / count_;
  }

  size_t count() const { return count_; }
  Field sum() const { return sum_; }
};

template <typename T>
class GeometricMean {
  T product_;
  size_t count_;

 public:
  GeometricMean() : product_(1), count_(0) {}

  void record(T value) {
    product_ *= value;
    ++count_;
  }

  std::optional<T> get() const {
    if (count_ == 0) {
      return std::nullopt;
    }

    return std::pow(product_, 1. / count_);
  }
  bool has_value() const { return count_ != 0; }

  T operator*() const {
    assert(count_ > 0);
    return std::pow(product_, 1. / count_);
  }

  size_t count() const { return count_; }
  T product() const { return product_; }
};

template <typename T, typename Comparator = std::less<T>>
class Minimum {
  std::optional<T> result_;

  [[no_unique_address]]
  Comparator comparator_;

 public:
  Minimum() : result_(std::nullopt) {}

  void record(T value) {
    if (!result_ || comparator_(value, *result_)) {
      result_ = value;
    }
  }

  void record(std::optional<T> value) {
    if (value) {
      record(*value);
    }
  }

  std::optional<T> get() const { return result_; }
  bool has_value() const { return result_.has_value(); }

  T operator*() const {
    assert(has_value());
    return *result_;
  }
};

template <typename T, typename Comparator = std::less<T>>
class Maximum {
  std::optional<T> result_;

  [[no_unique_address]]
  Comparator comparator_;

 public:
  Maximum() : result_(std::nullopt) {}

  void record(T value) {
    if (!result_ || comparator_(*result_, value)) {
      result_ = value;
    }
  }

  void record(std::optional<T> value) {
    if (value) {
      record(*value);
    }
  }

  std::optional<T> get() const { return result_; }
  bool has_value() const { return result_.has_value(); }

  T operator*() const {
    assert(has_value());
    return *result_;
  }
};

// Calculates minimum i, s.t. a_i = min_j a_j
template <typename T, typename Comparator = std::less<T>>
class ArgMinimum {
 public:
  struct Result {
    size_t index;
    T min;
  };

 private:
  std::optional<Result> result_;

  [[no_unique_address]]
  Comparator comparator_;

 public:
  ArgMinimum() : result_(std::nullopt) {}

  void record(size_t index, T value) {
    if (!result_ || comparator_(value, result_->min) ||
        !comparator_(result_->min, value) && index < result_->index) {
      result_ = Result{index, value};
    }
  }

  void record(size_t index, std::optional<T> value) {
    if (value) {
      record(index, *value);
    }
  }

  std::optional<Result> get() const { return result_; }
  bool has_value() const { return result_.has_value(); }

  Result operator*() const {
    assert(has_value());
    return *result_;
  }
  const Result* operator->() const {
    assert(has_value());
    return std::addressof(*result_);
  }
};

// Calculates minimum i, s.t. a_i = max_j a_j
template <typename T, typename Comparator = std::less<T>>
class ArgMaximum {
 public:
  struct Result {
    size_t index;
    T max;
  };

 private:
  std::optional<Result> result_;

  [[no_unique_address]]
  Comparator comparator_;

 public:
  ArgMaximum() : result_(std::nullopt) {}

  void record(size_t index, T value) {
    if (!result_ || comparator_(result_->max, value) ||
        !comparator_(value, result_->max) && index < result_->index) {
      result_ = Result{index, value};
    }
  }

  void record(size_t index, std::optional<T> value) {
    if (value) {
      record(index, *value);
    }
  }

  std::optional<Result> get() const { return result_; }
  bool has_value() const { return result_.has_value(); }

  Result operator*() const {
    assert(has_value());
    return *result_;
  }

  const Result* operator->() const {
    assert(has_value());
    return std::addressof(*result_);
  }
};

template <typename T>
class KahanSum {
  T compensation_ = 0;
  T sum_ = 0;

 public:
  void add(T value) {
    const T y = value - compensation_;
    const T x = sum_ + y;

    compensation_ = (x - sum_) - y;
    sum_ = x;
  }

  T sum() const { return sum_; }
};
