#pragma once

#include "Types.h"
#include "linear/problem/MILPProblem.h"

namespace mps {

template <typename Field>
class ProblemGenerator {
 public:
  static MILPProblem<Field> generate(const MPSParsingState<Field>& state) {

  }
};

}  // namespace mps
