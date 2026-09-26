#pragma once

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
    .feasibility = 1e-9,
    .pivot = 1e-7,
    .suspicious_pivot = 1e-4,
};

}  // namespace simplex
