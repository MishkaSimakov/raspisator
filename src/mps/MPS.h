#pragma once

#include "MPSParser.h"
#include "ProblemGenerator.h"
#include "Types.h"
#include "linear/problem/MILPProblem.h"

namespace mps {

template <typename Field>
class MPS {
 public:
  static MILPProblem<Field> read(std::istream& is, Format format) {
    auto state = MPSParser<Field>::parse(is, format);

    return ProblemGenerator<Field>::generate(state);
  }
};

}  // namespace mps
