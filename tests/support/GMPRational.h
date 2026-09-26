#pragma once

#include <gmpxx.h>

#include <cctype>
#include <charconv>
#include <compare>
#include <cstdlib>
#include <format>
#include <optional>
#include <ostream>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

#include "field/FieldTraits.h"
#include "simplex/Tolerance.h"

// Thin wrapper around mpq_class. mpq_class uses expression templates: `a - b`
// is a lazy __gmp_expr<...>, not an mpq_class. That breaks generic code which
// deduces types (`auto`, CTAD, ternaries, overload resolution on Field), so
// every operation here returns a concrete GMPRational instead.
class GMPRational {
  mpq_class value_;

 public:
  GMPRational() = default;

  template <typename T>
    requires std::is_arithmetic_v<T>
  GMPRational(T value) {
    // gmpxx has no constructors for long long, so pick a type it supports
    if constexpr (std::is_floating_point_v<T>) {
      value_ = static_cast<double>(value);
    } else if constexpr (std::is_signed_v<T>) {
      value_ = static_cast<long>(value);
    } else {
      value_ = static_cast<unsigned long>(value);
    }
  }

  explicit GMPRational(mpq_class value) : value_(std::move(value)) {
    value_.canonicalize();
  }

  const mpq_class& get() const { return value_; }

  GMPRational& operator+=(const GMPRational& other) {
    value_ += other.value_;
    return *this;
  }

  GMPRational& operator-=(const GMPRational& other) {
    value_ -= other.value_;
    return *this;
  }

  GMPRational& operator*=(const GMPRational& other) {
    value_ *= other.value_;
    return *this;
  }

  GMPRational& operator/=(const GMPRational& other) {
    value_ /= other.value_;
    return *this;
  }

  GMPRational operator-() const { return GMPRational(mpq_class(-value_)); }

  explicit operator double() const { return value_.get_d(); }

  friend GMPRational operator+(GMPRational first, const GMPRational& second) {
    return first += second;
  }

  friend GMPRational operator-(GMPRational first, const GMPRational& second) {
    return first -= second;
  }

  friend GMPRational operator*(GMPRational first, const GMPRational& second) {
    return first *= second;
  }

  friend GMPRational operator/(GMPRational first, const GMPRational& second) {
    return first /= second;
  }

  friend bool operator==(const GMPRational& first, const GMPRational& second) {
    return first.value_ == second.value_;
  }

  friend std::strong_ordering operator<=>(const GMPRational& first,
                                          const GMPRational& second) {
    return cmp(first.value_, second.value_) <=> 0;
  }

  friend GMPRational abs(const GMPRational& value) {
    mpq_class result;
    mpq_abs(result.get_mpq_t(), value.value_.get_mpq_t());

    return GMPRational(std::move(result));
  }

  friend std::string to_string(const GMPRational& value) {
    return value.value_.get_str();
  }

  friend std::ostream& operator<<(std::ostream& os, const GMPRational& value) {
    return os << value.value_;
  }
};

template <>
struct FieldTraits<GMPRational> {
  inline static const GMPRational tolerance = 0;

  static GMPRational abs(const GMPRational& value) {
    return is_strictly_negative(value) ? -value : value;
  }

  static GMPRational floor(const GMPRational& value) {
    mpz_class floored;
    mpz_fdiv_q(floored.get_mpz_t(), value.get().get_num_mpz_t(),
               value.get().get_den_mpz_t());

    return GMPRational(mpq_class(floored));
  }

  static GMPRational fractional(const GMPRational& value) {
    return value - floor(value);
  }

  static bool is_strictly_positive(const GMPRational& value) {
    return sgn(value.get()) > 0;
  }
  static bool is_strictly_negative(const GMPRational& value) {
    return sgn(value.get()) < 0;
  }

  static bool is_nonzero(const GMPRational& value) {
    return sgn(value.get()) != 0;
  }
  static bool should_drop(const GMPRational& value) {
    return sgn(value.get()) == 0;
  }

  static GMPRational exp2(int exponent) {
    mpq_class result = 1;

    if (exponent >= 0) {
      mpq_mul_2exp(result.get_mpq_t(), result.get_mpq_t(), exponent);
    } else {
      mpq_div_2exp(result.get_mpq_t(), result.get_mpq_t(), -exponent);
    }

    return GMPRational(std::move(result));
  }

  // Parses decimal numbers like "-12.5e-3" exactly.
  // Returns std::nullopt if parsing failed
  static std::optional<GMPRational> from_string(std::string_view string) {
    size_t pos = 0;

    bool negative = false;
    if (pos < string.size() && (string[pos] == '+' || string[pos] == '-')) {
      negative = string[pos] == '-';
      ++pos;
    }

    std::string digits;
    long exponent = 0;
    bool seen_point = false;

    for (; pos < string.size(); ++pos) {
      const char c = string[pos];

      if (std::isdigit(static_cast<unsigned char>(c)) != 0) {
        digits.push_back(c);

        if (seen_point) {
          --exponent;
        }
      } else if (c == '.' && !seen_point) {
        seen_point = true;
      } else {
        break;
      }
    }

    if (digits.empty()) {
      return std::nullopt;
    }

    if (pos < string.size() && (string[pos] == 'e' || string[pos] == 'E')) {
      ++pos;

      bool negative_exponent = false;
      if (pos < string.size() && (string[pos] == '+' || string[pos] == '-')) {
        negative_exponent = string[pos] == '-';
        ++pos;
      }

      unsigned long written_exponent;
      const char* begin = string.data() + pos;
      const char* end = string.data() + string.size();
      const auto [ptr, ec] = std::from_chars(begin, end, written_exponent);

      if (ec != std::errc{}) {
        return std::nullopt;
      }

      pos = ptr - string.data();
      exponent += negative_exponent ? -static_cast<long>(written_exponent)
                                    : static_cast<long>(written_exponent);
    }

    if (pos != string.size()) {
      return std::nullopt;
    }

    mpz_class power;
    mpz_ui_pow_ui(power.get_mpz_t(), 10, std::labs(exponent));

    mpq_class result(mpz_class(digits, 10));

    if (exponent >= 0) {
      result *= power;
    } else {
      result /= power;
    }

    if (negative) {
      result = -result;
    }

    return GMPRational(std::move(result));
  }
};

template <>
struct std::formatter<GMPRational, char> : std::formatter<std::string> {
  template <class FmtContext>
  auto format(const GMPRational& value, FmtContext& ctx) const {
    return std::formatter<std::string>::format(to_string(value), ctx);
  }
};

namespace simplex {

template <>
inline const Tolerance<GMPRational> kDefaultTolerance<GMPRational> = {
    .feasibility = 0,
    .pivot = 0,
    .suspicious_pivot = 0,
};

}  // namespace simplex
