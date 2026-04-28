#pragma once

#include <cmath>

template <typename Field>
struct FieldTraits;

template <>
struct FieldTraits<double> {
  constexpr static double tolerance = 1e-10;
  constexpr static double drop_tolerance = 1e-14;

  static double abs(double value) { return std::abs(value); }

  static double floor(double value) {
    double floored = std::floor(value);

    if (std::abs(value - floored - 1) < tolerance) {
      return floored + 1;
    }

    return floored;
  }

  static double fractional(double value) {
    double floored = std::floor(value);

    if (std::abs(value - floored - 1) < tolerance) {
      return 0;
    }

    return value - floored;
  }

  static bool is_strictly_negative(double value) { return value < -tolerance; }
  static bool is_strictly_positive(double value) { return value > tolerance; }

  static bool is_nonzero(double value) {
    return is_strictly_negative(value) || is_strictly_positive(value);
  }

  static bool should_drop(double value) {
    return std::abs(value) < drop_tolerance;
  }

  static double exp2(int exponent) { return std::exp2(exponent); }
};

// template <std::signed_integral T>
// struct FieldTraits<T> {
//   static T abs(T value) { return std::abs(value); }
//
//   static T floor(T value) { return value; }
//
//   static T fractional(T value) { return 0; }
//
//   static bool is_strictly_negative(T value) { return value < 0; }
//   static bool is_strictly_positive(T value) { return value > 0; }
//
//   static bool is_nonzero(T value) { return value != 0; }
// };
