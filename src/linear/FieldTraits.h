#pragma once

template <typename Field>
struct FieldTraits;

template <>
struct FieldTraits<double> {
  static double abs(double value) { return std::abs(value); }

  static double floor(double value) { return std::floor(value); }

  static double is_strictly_negative(double value) { return value < -1e-10; }
  static double is_strictly_positive(double value) { return value > 1e-10; }

  static bool is_nonzero(double value) {
    return is_strictly_negative(value) || is_strictly_positive(value);
  }
};
