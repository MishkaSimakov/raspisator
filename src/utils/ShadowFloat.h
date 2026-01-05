#pragma once
#include <format>
#include <type_traits>

#include "linear/BigInteger.h"
#include "linear/FieldTraits.h"

template <typename Field>
  requires std::is_floating_point_v<Field>
class ShadowFloat {
  Field value_;
  Rational shadow_;

  void precision_guard() {}

 public:
  ShadowFloat() = default;

  // NOLINTNEXTLINE(google-explicit-constructor)
  ShadowFloat(Field value) : value_(value), shadow_(value) {}

  ShadowFloat(Field value, Rational shadow)
      : value_(value), shadow_(std::move(shadow)) {}

  Rational get_shadow() const { return shadow_; }
  Field get_value() const { return value_; }

  operator Field() const { return value_; }

  // comparisons
  bool is_strictly_positive() const {
    bool value_result = FieldTraits<Field>::is_strictly_positive(value_);
    bool shadow_result = FieldTraits<Rational>::is_strictly_positive(shadow_);

    if (value_result != shadow_result) {
      throw std::runtime_error(
          "Comparison result is altered by rounding errors");
    }

    return value_result;
  }

  bool is_strictly_negative() const {
    bool value_result = FieldTraits<Field>::is_strictly_negative(value_);
    bool shadow_result = FieldTraits<Rational>::is_strictly_negative(shadow_);

    if (value_result != shadow_result) {
      throw std::runtime_error(
          "Comparison result is altered by rounding errors");
    }

    return value_result;
  }

  bool is_nonzero() const {
    bool value_result = FieldTraits<Field>::is_nonzero(value_);
    bool shadow_result = FieldTraits<Rational>::is_nonzero(shadow_);

    if (value_result != shadow_result) {
      throw std::runtime_error(
          "Comparison result is altered by rounding errors");
    }

    return value_result;
  }

  ShadowFloat& operator+=(const ShadowFloat& value) {
    value_ += value.get_value();
    shadow_ += value.get_shadow();

    precision_guard();

    return *this;
  }

  ShadowFloat& operator-=(const ShadowFloat& value) {
    value_ -= value.get_value();
    shadow_ -= value.get_shadow();

    precision_guard();

    return *this;
  }

  ShadowFloat& operator*=(const ShadowFloat& value) {
    value_ *= value.get_value();
    shadow_ *= value.get_shadow();

    precision_guard();

    return *this;
  }

  ShadowFloat& operator/=(const ShadowFloat& value) {
    value_ /= value.get_value();
    shadow_ /= value.get_shadow();

    precision_guard();

    return *this;
  }
};

template <typename Field>
ShadowFloat<Field> operator+(const ShadowFloat<Field>& left,
                             const ShadowFloat<Field>& right) {
  auto result = left;
  result += right;
  return result;
}

template <typename Field>
ShadowFloat<Field> operator-(const ShadowFloat<Field>& left,
                             const ShadowFloat<Field>& right) {
  auto result = left;
  result -= right;
  return result;
}

template <typename Field>
ShadowFloat<Field> operator*(const ShadowFloat<Field>& left,
                             const ShadowFloat<Field>& right) {
  auto result = left;
  result *= right;
  return result;
}

template <typename Field>
ShadowFloat<Field> operator/(const ShadowFloat<Field>& left,
                             const ShadowFloat<Field>& right) {
  auto result = left;
  result /= right;
  return result;
}

template <typename Field>
struct FieldTraits<ShadowFloat<Field>> {
  constexpr static double kEpsilon = 1e-10;

  static ShadowFloat<Field> abs(ShadowFloat<Field> value) {
    return {
        FieldTraits<Field>::abs(value.get_value()),
        FieldTraits<Rational>::abs(value.get_shadow()),
    };
  }

  static ShadowFloat<Field> floor(ShadowFloat<Field> value) {
    return {
        FieldTraits<Field>::floor(value.get_value()),
        FieldTraits<Rational>::floor(value.get_shadow()),
    };
  }

  static ShadowFloat<Field> fractional(ShadowFloat<Field> value) {
    return {
        FieldTraits<Field>::fractional(value.get_value()),
        FieldTraits<Rational>::fractional(value.get_shadow()),
    };
  }

  static bool is_strictly_negative(ShadowFloat<Field> value) {
    return value.is_strictly_negative();
  }
  static bool is_strictly_positive(ShadowFloat<Field> value) {
    return value.is_strictly_positive();
  }

  static bool is_nonzero(ShadowFloat<Field> value) {
    return value.is_nonzero();
  }
};

template <typename Field>
struct std::formatter<ShadowFloat<Field>, char> : std::formatter<std::string> {
  template <class FmtContext>
  auto format(const ShadowFloat<Field>& value, FmtContext& ctx) const {
    auto str = std::to_string(static_cast<Field>(value));

    return std::ranges::copy(str, ctx.out()).out;
  }
};
