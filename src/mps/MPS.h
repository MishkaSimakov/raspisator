#pragma once

#include "Format.h"
#include "detail/MPSParser.h"
#include "detail/ProblemGenerator.h"

#include "linear/problem/MILPProblem.h"

namespace mps {

template <typename Field>
MILPProblem<Field> read(std::istream& is, Format format) {
  auto state = detail::MPSParser<Field>::parse(is, format);

  return detail::ProblemGenerator<Field>::generate(state);
}

}  // namespace mps
