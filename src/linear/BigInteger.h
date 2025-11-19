#pragma once

#include <fmt/ostream.h>

#include <algorithm>
#include <compare>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "FieldTraits.h"

using std::vector, std::strong_ordering, std::string;

// TODO: Small Object Optimization
class BigInteger {
 private:
  bool is_positive_;

  vector<int> chunks_;

  static constexpr size_t kBase = 10;
  static constexpr size_t kChunkWidth = 9;
  static constexpr int kChunkMax = 1e9;

  static constexpr size_t kDifferenceForSimpleMultiplication = 8;

  void reverse(size_t shift) {
    if (chunks_.size() == 1) {
      return;
    }

    std::reverse(chunks_.begin(), chunks_.end());

    // shift chunks
    int modulo_pow = basePow(shift);
    int multiply_pow = basePow(kChunkWidth - shift);

    if (shift == 0) {
      return;
    }

    for (size_t i = 1; i < chunks_.size(); ++i) {
      chunks_[i - 1] += chunks_[i] % modulo_pow * multiply_pow;
      chunks_[i] /= modulo_pow;
    }
  }

  // compare absolute values
  static strong_ordering compareAbs(const BigInteger& first,
                                    const BigInteger& second) {
    return compareAbs(first.chunks_.data(), second.chunks_.data(), first.size(),
                      second.size());
  }

  template <class T>
  static strong_ordering compareAbs(const T* first, const T* second,
                                    size_t first_size, size_t second_size) {
    strong_ordering result = first_size <=> second_size;

    if (result != strong_ordering::equal) {
      return result;
    }

    // sizes are equal
    for (size_t i = 0; i < first_size; ++i) {
      result = first[first_size - i - 1] <=> second[first_size - i - 1];

      if (result != std::strong_ordering::equal) {
        return result;
      }
    }

    return strong_ordering::equal;
  }

  static void simpleMultiplication(const long long* first,
                                   const long long* second, long long* result,
                                   size_t count, int sign_multiplier) {
    for (size_t i = 0; i < count; ++i) {
      for (size_t j = 0; j < count; ++j) {
        result[i + j] += first[i] * second[j] * sign_multiplier;
      }
    }

    propagateCarry(result, count * 2);
  }

  // first, second - multiplicand and multiplier
  // result - array for result (size at least 2 * count)
  // count - size of first and second
  // sign_multiplier - if 1 result of multiplication will be added to result
  //                   if -1 - subtracted
  static void multiply(const long long* first, const long long* second,
                       long long* result, size_t count, int sign_multiplier) {
    if (count <= kDifferenceForSimpleMultiplication) {
      simpleMultiplication(first, second, result, count, sign_multiplier);
      return;
    }

    size_t half = count / 2;

    long long* temp_result = new long long[count]();
    // bd
    multiply(first, second, temp_result, half, sign_multiplier);

    for (size_t i = 0; i < half; ++i) {
      result[i] += temp_result[i];
      result[i + half] += temp_result[i + half] + temp_result[i];
      result[i + half * 2] += temp_result[i + half];
    }

    // ac
    std::fill(temp_result, temp_result + count, 0);
    multiply(first + half, second + half, temp_result, half, sign_multiplier);

    for (size_t i = 0; i < half; ++i) {
      result[i + half] += temp_result[i];
      result[i + half * 2] += temp_result[i] + temp_result[i + half];
      result[i + half * 3] += temp_result[i + half];
    }

    // (a - b)(d - c) sign
    int ab_sign_multiplier =
        compareAbs(first, first + half, half, half) == strong_ordering::greater
            ? 1
            : -1;
    int cd_sign_multiplier = compareAbs(second + half, second, half, half) ==
                                     strong_ordering::greater
                                 ? 1
                                 : -1;

    for (size_t i = 0; i < half; ++i) {
      temp_result[i] = (first[i] - first[i + half]) * ab_sign_multiplier;
      temp_result[i + half] =
          (second[i + half] - second[i]) * cd_sign_multiplier;
    }

    propagateCarry(temp_result, half);
    propagateCarry(temp_result + half, half);

    multiply(temp_result, temp_result + half, result + half, half,
             ab_sign_multiplier * cd_sign_multiplier * sign_multiplier);

    // propagateCarry result
    propagateCarry(result, count * 2);

    delete[] temp_result;
  }

  // calculate abs(a) - abs(b)
  void subtractAbs(const BigInteger& value) {
    size_t first_size = size();
    size_t second_size = value.size();

    chunks_.resize(std::max(first_size, second_size) + 1);

    const int* first = chunks_.data();
    const int* second = value.chunks_.data();

    if (compareAbs(first, second, first_size, second_size) ==
        std::strong_ordering::less) {
      std::swap(first, second);
      std::swap(first_size, second_size);

      is_positive_ = !is_positive_;
    }

    for (size_t i = 0; i < second_size; ++i) {
      chunks_[i] = first[i] - second[i];
    }
    for (size_t i = second_size; i < first_size; ++i) {
      chunks_[i] = first[i];
    }

    propagateCarry();

    // remove unnecessary chunks
    chunks_.resize(sizeWithoutLeadingZeros());

    // prevent -0
    if (!*this) {
      is_positive_ = true;
    }
  }

  void propagateCarry() {
    propagateCarry(chunks_.data(), size());

    int back = chunks_.back();
    if (back >= kChunkMax) {
      chunks_.back() = back % kChunkMax;
      chunks_.push_back(back / kChunkMax);
    }
  }

  template <class T>
  static void propagateCarry(T* values, size_t count) {
    long long carry = 0;

    for (size_t i = 0; i < count; ++i) {
      carry += static_cast<long long>(values[i]);

      if (i == count - 1) {
        values[i] = carry;
        break;
      }

      if (carry < 0) {
        long long new_carry = (carry + 1) / kChunkMax - 1;

        values[i] = static_cast<T>(carry + kChunkMax * (-new_carry));

        carry = new_carry;
      } else {
        values[i] = static_cast<T>(carry % kChunkMax);
        carry /= kChunkMax;
      }
    }
  }

  size_t sizeWithoutLeadingZeros() const {
    return sizeWithoutLeadingZeros(chunks_.data(), size());
  }

  template <class T>
  static size_t sizeWithoutLeadingZeros(const T* values, size_t count) {
    while (count > 1 && values[count - 1] == 0) {
      --count;
    }

    return count;
  }

  // calculate abs(a) + abs(b)
  void addAbs(const BigInteger& value) {
    size_t addend_size = value.size();

    chunks_.resize(std::max(size(), addend_size), 0);

    for (size_t i = 0; i < addend_size; ++i) {
      chunks_[i] += value.chunks_[i];
    }

    propagateCarry();
  }

  size_t size() const { return chunks_.size(); }

  static size_t ceilToPow2(size_t value) {
    size_t result = 1;
    while (result < value) {
      result *= 2;
    }

    return result;
  }

  // methods for division
  struct DivisionDTO {
    long long* dividend;
    long long* divisor;
    size_t dividend_size;
    size_t divisor_size;

    long long multiplier;
  };

  static bool isNegative(const long long* value, size_t size) {
    size = sizeWithoutLeadingZeros(value, size);
    return value[size - 1] < 0;
  }

  void divideAbs(DivisionDTO& division_dto, bool save_quotient) {
    long long* dividend = division_dto.dividend;
    long long* divisor = division_dto.divisor;

    size_t& dividend_size = division_dto.dividend_size;
    size_t& divisor_size = division_dto.divisor_size;

    if (save_quotient) {
      chunks_.clear();
      chunks_.resize(dividend_size - divisor_size + 1, 0);
    }

    size_t divisor_shift = dividend_size - divisor_size + 1;

    do {
      divisor_shift = std::min(divisor_shift - 1, dividend_size - divisor_size);

      if (dividend_size == divisor_size + divisor_shift &&
          compareAbs(dividend + divisor_shift, divisor,
                     dividend_size - divisor_shift,
                     divisor_size) != strong_ordering::less) {
        for (size_t i = 0; i < divisor_size; ++i) {
          dividend[i + divisor_shift] -= divisor[i];
        }

        propagateCarry(dividend + divisor_shift, dividend_size - divisor_shift);

        dividend_size = sizeWithoutLeadingZeros(dividend, dividend_size);

        if (save_quotient) {
          chunks_[divisor_shift] += 1;
        }

        continue;
      }

      if (dividend_size == divisor_size + divisor_shift) {
        continue;
      }

      // something like Knuth algorithm for division
      // approximation <= base^2 / (base/2) = 2 * base
      long long approximation = (dividend[dividend_size - 1] * kChunkMax +
                                 dividend[dividend_size - 2]) /
                                divisor[divisor_size - 1];

      // calculate dividend - divisor * approximation
      for (size_t i = 0; i < divisor_size; ++i) {
        dividend[i + divisor_shift] -= divisor[i] * approximation;
      }

      propagateCarry(dividend + divisor_shift, dividend_size - divisor_shift);

      while (isNegative(dividend, dividend_size)) {
        approximation -= 1;

        for (size_t i = 0; i < divisor_size; ++i) {
          dividend[i + divisor_shift] += divisor[i];
        }

        propagateCarry(dividend + divisor_shift, dividend_size - divisor_shift);
      }

      dividend_size = sizeWithoutLeadingZeros(dividend, dividend_size);

      if (save_quotient) {
        chunks_[divisor_shift] += static_cast<int>(approximation % kChunkMax);
      }
    } while (divisor_shift != 0 && dividend_size >= divisor_size);

    if (save_quotient) {
      chunks_.resize(sizeWithoutLeadingZeros());
    }
  }

  static DivisionDTO getDivisionDTO(const BigInteger& dividend,
                                    const BigInteger& divisor) {
    size_t dividend_size = dividend.size();
    size_t divisor_size = divisor.size();

    long long first_digit = divisor.chunks_.back();

    long long multiplier =
        (first_digit * 2 < kChunkMax) ? kChunkMax / (first_digit + 1) : 1;

    long long* dividend_ptr = new long long[dividend_size + 1];
    long long* divisor_ptr = new long long[divisor_size];

    for (size_t i = 0; i < dividend_size; ++i) {
      dividend_ptr[i] = multiplier * dividend.chunks_[i];
    }
    dividend_ptr[dividend_size] = 0;

    for (size_t i = 0; i < divisor_size; ++i) {
      divisor_ptr[i] = multiplier * divisor.chunks_[i];
    }

    propagateCarry(dividend_ptr, dividend_size + 1);
    propagateCarry(divisor_ptr, divisor_size);

    if (dividend_ptr[dividend_size] != 0) {
      ++dividend_size;
    }

    return {dividend_ptr, divisor_ptr, dividend_size, divisor_size, multiplier};
  }

 public:
  BigInteger(unsigned long long value, bool is_positive)
      : is_positive_(is_positive || value == 0) {
    if (value == 0) {
      chunks_.push_back(0);
      return;
    }

    while (value != 0) {
      chunks_.push_back(static_cast<int>(value % kChunkMax));
      value /= kChunkMax;
    }
  }

  BigInteger(long long value = 0) : BigInteger(std::abs(value), value >= 0) {}

  BigInteger(const char* string, unsigned long long length) {
    chunks_.reserve(length / kChunkWidth);

    std::stringstream ss(string);
    ss >> *this;
  }

  explicit operator bool() const {
    return chunks_.size() != 1 || chunks_.front() != 0;
  }

  BigInteger& operator+=(const BigInteger& value) {
    if (value.is_positive_ != is_positive_) {
      // negative += positive
      // or positive += negative
      subtractAbs(value);
    } else {
      // positive += positive
      // or negative += negative
      addAbs(value);
    }

    return *this;
  }

  BigInteger& operator-=(const BigInteger& value) {
    if (value.is_positive_ != is_positive_) {
      // negative -= positive
      // or positive -= negative
      addAbs(value);
    } else {
      // positive -= positive
      // or negative -= negative
      subtractAbs(value);
    }

    return *this;
  }

  BigInteger& operator*=(const BigInteger& value) {
    if (!*this || !value) {
      zero();

      return *this;
    }

    size_t new_size = ceilToPow2(std::max(size(), value.size()));

    vector<long long> first_copy(new_size);
    vector<long long> second_copy(new_size);
    vector<long long> result(new_size * 2);

    std::copy(chunks_.begin(), chunks_.end(), first_copy.begin());
    std::copy(value.chunks_.begin(), value.chunks_.end(), second_copy.begin());

    BigInteger::multiply(first_copy.data(), second_copy.data(), result.data(),
                         new_size, 1);

    size_t real_size = sizeWithoutLeadingZeros(result.data(), result.size());

    chunks_.resize(real_size);
    std::copy(result.begin(), result.begin() + real_size, chunks_.begin());

    is_positive_ = is_positive_ == value.is_positive_;

    return *this;
  }

  BigInteger& operator/=(const BigInteger& value) {
    if (size() < value.size()) {
      zero();
      return *this;
    }

    DivisionDTO division_dto = getDivisionDTO(*this, value);
    divideAbs(division_dto, true);

    delete[] division_dto.divisor;
    delete[] division_dto.dividend;

    is_positive_ = !(*this) || is_positive_ == value.is_positive_;

    return *this;
  }

  BigInteger& operator%=(const BigInteger& modulo) {
    if (size() < modulo.size()) {
      return *this;
    }

    DivisionDTO division_dto = getDivisionDTO(*this, modulo);
    divideAbs(division_dto, true);

    chunks_.resize(division_dto.dividend_size);

    std::copy(division_dto.dividend,
              division_dto.dividend + division_dto.dividend_size,
              chunks_.begin());

    *this /= division_dto.multiplier;

    delete[] division_dto.divisor;
    delete[] division_dto.dividend;

    return *this;
  }

  BigInteger& operator++() { return *this += 1; }

  BigInteger& operator--() { return *this -= 1; }

  BigInteger operator++(int) {
    BigInteger copy(*this);

    *this += 1;

    return copy;
  }

  BigInteger operator--(int) {
    BigInteger copy(*this);

    *this -= 1;

    return copy;
  }

  friend std::ostream& operator<<(std::ostream& os, const BigInteger& value);

  friend std::istream& operator>>(std::istream& is, BigInteger& value);

  BigInteger operator-() const {
    BigInteger copy(*this);
    copy.negate();

    return copy;
  }

  std::strong_ordering operator<=>(const BigInteger& other) const {
    if (!is_positive_ && other.is_positive_) {
      return std::strong_ordering::less;
    }
    if (is_positive_ && !other.is_positive_) {
      return std::strong_ordering::greater;
    }

    strong_ordering result = compareAbs(*this, other);

    if (result == strong_ordering::equal) {
      return strong_ordering::equal;
    }

    if (!is_positive_) {
      return result == strong_ordering::less ? strong_ordering::greater
                                             : strong_ordering::less;
    }

    return result;
  }

  bool operator==(const BigInteger& other) const = default;

  void zero() {
    is_positive_ = true;
    chunks_.clear();
    chunks_.push_back(0);
  }

  string toString() const {
    std::stringstream stream;
    stream << *this;

    return stream.str();
  }

  void negate() {
    if (*this) {
      is_positive_ = !is_positive_;
    }
  }

  // shifts number by {offset} digits
  BigInteger& operator<<=(size_t offset) {
    // zero doesn't need to be shifted
    if (!*this) {
      return *this;
    }

    size_t chunks_offset = offset / kChunkWidth;
    size_t old_size = size();

    chunks_.resize(old_size + chunks_offset);

    for (size_t i = 0; i < old_size; ++i) {
      size_t index = old_size - i - 1;
      chunks_[index + chunks_offset] = chunks_[index];
    }
    for (size_t i = 0; i < chunks_offset; ++i) {
      chunks_[i] = 0;
    }

    offset %= kChunkWidth;

    if (offset == 0) {
      return *this;
    }

    old_size = size();
    int multiplier = basePow(offset);
    int modulo = basePow(kChunkWidth - offset);

    if (chunks_.back() >= modulo) {
      chunks_.push_back(0);
    }

    for (size_t i = 0; i < old_size; ++i) {
      size_t index = old_size - i - 1;
      if (chunks_[index] >= modulo) {
        chunks_[index + 1] += chunks_[index] / modulo;
      }

      chunks_[index] %= modulo;
      chunks_[index] *= multiplier;
    }

    return *this;
  }

  static int basePow(size_t power) {
    return static_cast<int>(std::pow(kBase, power));
  }

  static int baseLog(int value) { return static_cast<int>(std::log10(value)); }
};

inline BigInteger operator+(const BigInteger& first, const BigInteger& second) {
  BigInteger copy(first);

  copy += second;

  return copy;
}

inline BigInteger operator-(const BigInteger& first, const BigInteger& second) {
  BigInteger copy(first);

  copy -= second;

  return copy;
}

inline BigInteger operator*(const BigInteger& first, const BigInteger& second) {
  BigInteger copy(first);

  copy *= second;

  return copy;
}

inline BigInteger operator/(const BigInteger& first, const BigInteger& second) {
  BigInteger copy(first);

  copy /= second;

  return copy;
}

inline BigInteger operator%(const BigInteger& first, const BigInteger& second) {
  BigInteger copy(first);

  copy %= second;

  return copy;
}

inline std::ostream& operator<<(std::ostream& os, const BigInteger& value) {
  if (!value.is_positive_) {
    os << "-";
  }

  size_t chunks_count = value.size();

  for (size_t i = 0; i < chunks_count; ++i) {
    if (i != 0) {
      os << std::setfill('0') << std::setw(BigInteger::kChunkWidth);
    }

    os << value.chunks_[chunks_count - i - 1];
  }

  return os;
}

inline std::istream& operator>>(std::istream& is, BigInteger& value) {
  vector<int>& chunks = value.chunks_;

  chunks.clear();
  value.is_positive_ = true;

  // skip whitespaces at begin
  while (std::isspace(is.peek()) != 0) {
    is.get();
  }

  // read sign
  int first_char = is.peek();
  if (first_char == '-' || first_char == '+') {
    value.is_positive_ = first_char == '+';
    is.get();
  }

  int chunk = 0;
  int position_in_chunk = BigInteger::kChunkMax;

  bool non_zero = false;

  while ('0' <= is.peek() && is.peek() <= '9') {
    int input = is.get();

    if (input != '0') {
      non_zero = true;
    }

    position_in_chunk /= BigInteger::kBase;

    chunk += (input - '0') * position_in_chunk;

    if (position_in_chunk == 1) {
      chunks.push_back(chunk);
      chunk = 0;
      position_in_chunk = BigInteger::kChunkMax;
    }
  }

  if (position_in_chunk != BigInteger::kChunkMax) {
    chunk /= position_in_chunk;
    chunks.push_back(chunk);
  }

  // if stream was empty fill with 0
  if (chunks.empty()) {
    chunks.push_back(0);
  }

  // prevent -0
  if (!non_zero) {
    value.is_positive_ = true;
  }

  value.reverse(BigInteger::baseLog(position_in_chunk) %
                BigInteger::kChunkWidth);

  chunks.resize(value.sizeWithoutLeadingZeros());

  return is;
}

inline BigInteger operator""_bi(const char* string, unsigned long length) {
  return {string, length};
}

inline BigInteger operator""_bi(unsigned long long value) {
  return {value, true};
}

// Rational
class Rational {
 private:
  static constexpr size_t kPrecisionForDouble = 25;

  BigInteger numerator_;
  BigInteger denominator_;

  static BigInteger gcd(BigInteger first, BigInteger second) {
    if (first < second) {
      std::swap(first, second);
    }

    while (second) {
      first %= second;
      std::swap(first, second);
    }

    return first;
  }

  void normalize() {
    // only one rational for zero - 0 / 1
    if (!numerator_) {
      denominator_ = 1;
      return;
    }
    // minus only in numerator
    if (denominator_ < 0) {
      denominator_.negate();
      numerator_.negate();
    }

    bool is_negative = false;
    if (numerator_ < 0) {
      is_negative = true;
      numerator_.negate();
    }

    BigInteger common_divider = Rational::gcd(numerator_, denominator_);

    numerator_ /= common_divider;
    denominator_ /= common_divider;

    if (is_negative) {
      numerator_.negate();
    }
  }

 public:
  Rational(const BigInteger& value) : numerator_(value), denominator_(1) {}

  Rational(long long value = 0) : numerator_(value), denominator_(1) {}

  Rational floor() const { return numerator_ / denominator_; }

  Rational& operator+=(const Rational& other) {
    numerator_ *= other.denominator_;
    numerator_ += other.numerator_ * denominator_;

    denominator_ *= other.denominator_;

    normalize();

    return *this;
  }

  Rational& operator-=(const Rational& other) {
    numerator_ *= other.denominator_;
    numerator_ -= other.numerator_ * denominator_;
    denominator_ *= other.denominator_;

    normalize();

    return *this;
  }

  Rational& operator*=(const Rational& other) {
    numerator_ *= other.numerator_;
    denominator_ *= other.denominator_;

    normalize();

    return *this;
  }

  Rational& operator/=(const Rational& other) {
    numerator_ *= other.denominator_;
    denominator_ *= other.numerator_;

    normalize();

    return *this;
  }

  Rational operator-() const {
    Rational copy(*this);

    copy.numerator_.negate();

    return copy;
  }

  string toString() const {
    string result;

    result += numerator_.toString();

    if (denominator_ != 1) {
      result += '/';
      result += denominator_.toString();
    }

    return result;
  }

  strong_ordering operator<=>(const Rational& other) const {
    return numerator_ * other.denominator_ <=> other.numerator_ * denominator_;
  }

  bool operator==(const Rational& other) const = default;

  string asDecimal(size_t precision = 0) const {
    BigInteger numerator_copy(numerator_);

    bool negative = false;
    if (numerator_copy < 0) {
      negative = true;
      numerator_copy.negate();
    }

    numerator_copy <<= precision;

    numerator_copy /= denominator_;

    std::string result = numerator_copy.toString();

    if (result.size() < precision + 1) {
      result.insert(0, precision + 1 - result.size(), '0');
    }

    if (precision > 0) {
      result.insert(result.length() - precision, ".");
    }
    if (negative) {
      result.insert(0, "-");
    }

    return result;
  }

  explicit operator double() const {
    std::stringstream ss(asDecimal(kPrecisionForDouble));

    double result;
    ss >> result;

    return result;
  }
};

inline Rational operator+(const Rational& first, const Rational& second) {
  Rational copy(first);

  copy += second;

  return copy;
}

inline Rational operator-(const Rational& first, const Rational& second) {
  Rational copy(first);

  copy -= second;

  return copy;
}

inline Rational operator*(const Rational& first, const Rational& second) {
  Rational copy(first);

  copy *= second;

  return copy;
}

inline Rational operator/(const Rational& first, const Rational& second) {
  Rational copy(first);

  copy /= second;

  return copy;
}

inline std::ostream& operator<<(std::ostream& os, const Rational& rational) {
  os << rational.toString();

  return os;
}

inline std::istream& operator>>(std::istream& is, Rational& rational) {
  while (std::isspace(is.peek()) != 0) {
    is.get();
  }

  // read sign
  bool negate = false;
  int first_char = is.peek();
  if (first_char == '-' || first_char == '+') {
    negate = true;
    is.get();
  }

  rational = 0;
  bool fractional_part = false;
  size_t fractional_size = 0;

  while ('0' <= is.peek() && is.peek() <= '9' ||
         !fractional_part && is.peek() == '.') {
    int input = is.get();

    if (fractional_part) {
      ++fractional_size;
    }

    if (input == '.') {
      fractional_part = true;
    } else {
      rational *= 10;
      rational += input - '0';
    }
  }

  if (negate) {
    rational *= -1;
  }

  for (size_t i = 0; i < fractional_size; ++i) {
    rational /= 10;
  }

  return is;
}

template <>
struct FieldTraits<Rational> {
  static Rational floor(const Rational& value) { return value.floor(); }
};

template <>
struct std::formatter<Rational, char> {
  template <class ParseContext>
  constexpr ParseContext::iterator parse(ParseContext& ctx) {
    // auto it = ctx.begin();
    // if (it == ctx.end()) return it;
    //
    // if (*it == '#') {
    //   quoted = true;
    //   ++it;
    // }
    // if (it != ctx.end() && *it != '}')
    //   throw std::format_error("Invalid format args for QuotableString.");
    //
    return ctx.begin();
  }

  template <class FmtContext>
  FmtContext::iterator format(const Rational& value, FmtContext& ctx) const {
    std::ostringstream out;
    out << value;

    return std::ranges::copy(std::move(out).str(), ctx.out()).out;
  }
};
