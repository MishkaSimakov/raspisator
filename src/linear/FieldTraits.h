#pragma once

#include <cmath>

template <typename Field>
struct FieldTraits;

template <>
struct FieldTraits<double> {
  constexpr static double kEpsilon = 1e-10;

  static double abs(double value) { return std::abs(value); }

  static double floor(double value) {
    double floored = std::floor(value);

    if (std::abs(value - floored - 1) < kEpsilon) {
      return floored + 1;
    }

    return floored;
  }

  static double fractional(double value) {
    double floored = std::floor(value);

    if (std::abs(value - floored - 1) < kEpsilon) {
      return 0;
    }

    return value - floored;
  }

  static bool is_strictly_negative(double value) { return value < -kEpsilon; }
  static bool is_strictly_positive(double value) { return value > kEpsilon; }

  static bool is_nonzero(double value) {
    return is_strictly_negative(value) || is_strictly_positive(value);
  }
};

template <std::integral T>
struct FieldTraits<T> {
  static T abs(T value) { return std::abs(value); }

  static T floor(T value) { return value; }

  static T fractional(T value) { return 0; }

  static bool is_strictly_negative(T value) { return value < 0; }
  static bool is_strictly_positive(T value) { return value > 0; }

  static bool is_nonzero(T value) { return value != 0; }
};
