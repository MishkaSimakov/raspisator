#pragma once

#include "Format.h"
#include "detail/MPSParser.h"
#include "detail/ProblemGenerator.h"

#include "problem/MILP.h"

namespace mps {

template <typename Field>
problem::MILP<Field> read(std::istream& is, Format format) {
  auto state = detail::MPSParser<Field>::parse(is, format);

  return detail::ProblemGenerator<Field>::generate(state);
}

}  // namespace mps
