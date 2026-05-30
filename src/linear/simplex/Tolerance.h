#pragma once

#include "linear/BigInteger.h"

namespace simplex {

template <typename Field>
struct Tolerance {
  Field feasibility;
  Field pivot;
};

// Default tolerance values for standard fields
template <typename Field>
constexpr Tolerance<Field> kDefaultTolerance;

template <>
constexpr Tolerance kDefaultTolerance<double> = {
    .feasibility = 1e-7,
    .pivot = 1e-7,
};

template <>
constexpr Tolerance kDefaultTolerance<Rational> = {
    .feasibility = 0,
    .pivot = 0,
};

}  // namespace simplex
