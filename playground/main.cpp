#include <iostream>

#include "linear/BigInteger.h"
#include "linear/BranchAndBound.h"
#include "linear/MPS.h"
#include "linear/SimplexMethod.h"

using Field = Rational;

int main() {
  Matrix<Field> A = {{8000, 4000, 1, 0}, {15, 30, 0, 1}};

  Matrix<Field> b = {{40'000}, {200}};
  Matrix<Field> c = {{100, 150, 0, 0}};

  std::vector<size_t> integer = {0, 1};

  MILPProblem<Field> problem(A, b, c, integer);

  BranchAndBound<Field, SimplexMethod<Field>> solver(problem);
  auto [value, point] = solver.solve();

  std::cout << GraphvizBuilder<Field>().build(solver.get_root()) << std::endl;

  std::cout << point << std::endl;

  return 0;
}
