// Degeneracy in the simplex method (as well as in real life) is very dangerous.
// https://ocw.mit.edu/courses/15-053-optimization-methods-in-management-science-spring-2013/1542510abf20fa6145dfb0b1551be6c2_MIT15_053S13_tut07.pdf
// Therefore the implementation must be rigorously tested on problems with
// degeneracy.

#include <gtest/gtest.h>

#include "Assertions.h"
#include "linear/BigInteger.h"
#include "linear/Matrix.h"
#include "linear/SimplexMethod.h"

TEST(DegenerateSimplexMethodTests, SimpleTests) {
  Matrix<Rational> A = {
      {0, 1, 0},
      {-1, 0, 0},
  };

  Matrix<Rational> b = {{1}, {0}};
  Matrix<Rational> c = {{1, 0, 0}};

  auto bfs = SimplexMethod(A, b, c).find_bfs();

  ASSERT_TRUE(bfs.has_value());
  ASSERT_SETS_EQ(bfs->basic_variables, (std::vector<size_t>{0, 1}));
}
