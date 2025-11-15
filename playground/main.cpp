#include <iostream>

#include "linear/BigInteger.h"
#include "linear/MPS.h"
#include "linear/SimplexMethod.h"

int main() {
  auto reader = MPSReader<Rational>();
  reader.read("resources/adlittle.mps");

  auto [A, b, c] = reader.get_canonical_representation();

  std::cout << "A" << std::endl;
  std::cout << A << std::endl;

  std::cout << "b" << std::endl;
  std::cout << b << std::endl;

  std::cout << "c" << std::endl;
  std::cout << c << std::endl;

  for (size_t j = 0; j < A.get_width(); ++j) {
    bool all_zeros = true;

    for (size_t i = 0; i < A.get_height(); ++i) {
      if (A[i, j] != 0) {
        all_zeros = false;
        break;
      }
    }

    if (all_zeros) {
      std::cout << "there is zero column in A!" << std::endl;
    }
  }

  auto solver = SimplexMethod<Rational>(A, b, c);
  auto solution = solver.solve();

  if (std::holds_alternative<FiniteSolution<Rational>>(solution)) {
    std::cout << std::get<FiniteSolution<Rational>>(solution).point
              << std::endl;
    std::cout << std::get<FiniteSolution<Rational>>(solution).value
              << std::endl;
  } else if (std::holds_alternative<InfiniteSolution>(solution)) {
    std::cout << "Infinite solution" << std::endl;
  } else {
    std::cout << "No feasible points" << std::endl;
  }

  return 0;
}
