#pragma once

#include <optional>

#include "linear/FieldTraits.h"

enum class BoundViolationType {
  VIOLATE_UPPER_BOUND,
  VIOLATE_LOWER_BOUND,
  NONE
};

template <typename Field>
struct BoundViolation {
  BoundViolationType type;
  Field value;

  auto operator<=>(const BoundViolation& other) const {
    return value <=> other.value;
  }
};

template <typename Field>
struct Bound {
  std::optional<Field> lower = std::nullopt;
  std::optional<Field> upper = std::nullopt;

  Bound() = default;

  Bound(std::optional<Field> lower, std::optional<Field> upper)
      : lower(lower), upper(upper) {}

  Bound& operator+=(const Bound& other) {
    lower = lower && other.lower ? std::optional(*lower + *other.lower)
                                 : std::nullopt;
    upper = upper && other.upper ? std::optional(*upper + *other.upper)
                                 : std::nullopt;

    return *this;
  }

  Bound& operator-=(const Bound& other) {
    lower = lower && other.lower ? std::optional(*lower - *other.upper)
                                 : std::nullopt;
    upper = upper && other.upper ? std::optional(*upper - *other.lower)
                                 : std::nullopt;

    return *this;
  }

  Bound& operator*=(Field scalar) {
    if (scalar >= 0) {
      lower = lower.transform([scalar](Field value) { return value * scalar; });
      upper = upper.transform([scalar](Field value) { return value * scalar; });
    } else {
      auto old_lower = lower;
      lower = upper.transform([scalar](Field value) { return value * scalar; });
      upper =
          old_lower.transform([scalar](Field value) { return value * scalar; });
    }

    return *this;
  }

  Bound& operator/=(Field scalar) {
    *this *= Field(1) / scalar;

    return *this;
  }

  Bound& operator^=(const Bound& other) {
    stricten_lower(other.lower);
    stricten_upper(other.upper);

    return *this;
  }

  Bound operator-() const { return *this * Field(-1); }

  void stricten_lower(std::optional<Field> new_bound) {
    if (!new_bound) {
      return;
    }

    if (!lower) {
      lower = new_bound;
    } else {
      lower = std::max(*lower, *new_bound);
    }
  }

  void stricten_upper(std::optional<Field> new_bound) {
    if (!new_bound) {
      return;
    }

    if (!upper) {
      upper = new_bound;
    } else {
      upper = std::min(*upper, *new_bound);
    }
  }

  bool is_fixed() const {
    return lower && upper && !FieldTraits<Field>::is_nonzero(*lower - *upper);
  }

  bool is_infeasible() const {
    return lower && upper &&
           FieldTraits<Field>::is_strictly_negative(*upper - *lower);
  }

  // Returns a point inside bounds
  Field get_point_inside() const {
    if (lower) {
      return *lower;
    }

    if (upper) {
      return *upper;
    }

    return 0;
  }

  bool is_inside(Field value) const {
    if (lower && FieldTraits<Field>::is_strictly_negative(value - *lower)) {
      return false;
    }

    if (upper && FieldTraits<Field>::is_strictly_negative(*upper - value)) {
      return false;
    }

    return true;
  }

  BoundViolation<Field> get_violation(Field value) const {
    if (lower) {
      Field lower_violation = *lower - value;

      if (FieldTraits<Field>::is_strictly_positive(lower_violation)) {
        return {
            .type = BoundViolationType::VIOLATE_LOWER_BOUND,
            .value = lower_violation,
        };
      }
    }

    if (upper) {
      Field upper_violation = value - *upper;

      if (FieldTraits<Field>::is_strictly_positive(upper_violation)) {
        return {
            .type = BoundViolationType::VIOLATE_UPPER_BOUND,
            .value = upper_violation,
        };
      }
    }

    return {
        .type = BoundViolationType::NONE,
        .value = 0,
    };
  }
};

template <typename Field>
Bound<Field> operator+(const Bound<Field>& left, const Bound<Field>& right) {
  auto copy = left;
  copy += right;
  return copy;
}

template <typename Field>
Bound<Field> operator-(const Bound<Field>& left, const Bound<Field>& right) {
  auto copy = left;
  copy -= right;
  return copy;
}

template <typename Field>
Bound<Field> operator*(Field left, const Bound<Field>& right) {
  auto copy = right;
  copy *= left;
  return copy;
}

template <typename Field>
Bound<Field> operator*(const Bound<Field>& left, Field right) {
  return right * left;
}

template <typename Field>
Bound<Field> operator/(const Bound<Field>& left, Field right) {
  auto copy = left;
  copy /= right;
  return copy;
}

template <typename Field>
Bound<Field> operator^(const Bound<Field>& left, const Bound<Field>& right) {
  auto copy = left;
  copy ^= right;
  return copy;
}

template <typename Field>
std::ostream& operator<<(std::ostream& os, const Bound<Field>& bound) {
  os << "(";

  if (bound.lower) {
    os << *bound.lower;
  } else {
    os << "-∞";
  }

  os << ", ";

  if (bound.upper) {
    os << *bound.upper;
  } else {
    os << "+∞";
  }

  os << ")";

  return os;
}

template <typename Field>
struct std::formatter<Bound<Field>, char> : std::formatter<std::string> {
  template <class FmtContext>
  auto format(const Bound<Field>& value, FmtContext& ctx) const {
    std::ostringstream out;
    out << value;
    return std::ranges::copy(std::move(out).str(), ctx.out()).out;
  }
};
