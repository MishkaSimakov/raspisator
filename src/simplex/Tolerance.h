#pragma once

#include "field/BigInteger.h"

namespace simplex {

template <typename Field>
struct Tolerance {
  Field feasibility;
  Field pivot;
  Field suspicious_pivot;
};

// Default tolerance values for standard fields
template <typename Field>
const Tolerance<Field> kDefaultTolerance;

template <>
inline const Tolerance<double> kDefaultTolerance<double> = {
    .feasibility = 1e-7,
    .pivot = 1e-7,
    .suspicious_pivot = 1e-4,
};

template <>
inline const Tolerance<Rational> kDefaultTolerance<Rational> = {
    .feasibility = 0,
    .pivot = 0,
    .suspicious_pivot = 0,
};

}  // namespace simplex
