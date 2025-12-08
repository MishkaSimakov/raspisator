#pragma once

#include <cmath>

template <typename Field>
struct FieldTraits;

template <>
struct FieldTraits<double> {
  static double abs(double value) { return std::abs(value); }

  static double floor(double value) {
    double floored = std::floor(value);

    if (std::abs(value - floored - 1) < 1e-10) {
      return floored + 1;
    }

    return floored;
  }

  static double fractional(double value) {
    double floored = std::floor(value);

    if (std::abs(value - floored - 1) < 1e-10) {
      return 0;
    }

    return value - floored;
  }

  static bool is_strictly_negative(double value) { return value < -1e-10; }
  static bool is_strictly_positive(double value) { return value > 1e-10; }

  static bool is_nonzero(double value) {
    return is_strictly_negative(value) || is_strictly_positive(value);
  }
};
