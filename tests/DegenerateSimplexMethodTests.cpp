// Degeneracy in the simplex method (as well as in real life) is very dangerous.
// https://ocw.mit.edu/courses/15-053-optimization-methods-in-management-science-spring-2013/1542510abf20fa6145dfb0b1551be6c2_MIT15_053S13_tut07.pdf
// Therefore the implementation must be rigorously tested on problems with
// degeneracy.

#include <gtest/gtest.h>

#include "Assertions.h"
#include "linear/BigInteger.h"
#include "linear/BoundedSimplexMethod.h"
#include "linear/matrix/Matrix.h"
//
// TEST(DegenerateSimplexMethodTests, SimpleTest) {
//   Matrix<Rational> A = {
//       {0, 1, 0},
//       {-1, 0, 0},
//   };
//
//   Matrix<Rational> b = {{1}, {0}};
//   Matrix<Rational> c = {{1, 0, 0}};
//
//   auto bfs = BoundedSimplexMethod(CSCMatrix(A), b, c).find_bfs();
//
//   ASSERT_TRUE(bfs.has_value());
//   ASSERT_SETS_EQ(bfs->basic_variables, (std::vector<size_t>{0, 1}));
// }
//
// TEST(DegenerateSimplexMethodTests, ZerosInb1) {
//   size_t N = 20;
//
//   Matrix<Rational> A(N, N + 1, 0);
//
//   for (size_t i = 0; i < N; ++i) {
//     A[i, i] = -Rational(i + 1);
//   }
//
//   Matrix<Rational> b(N, 1, 0);
//   Matrix<Rational> c(1, N + 1, 1);
//
//   auto bfs = BoundedSimplexMethod(CSCMatrix(A), b, c).find_bfs();
//
//   ASSERT_TRUE(bfs.has_value());
//
//   std::vector<size_t> expected(N);
//   std::iota(expected.begin(), expected.end(), 0);
//
//   ASSERT_SETS_EQ(bfs->basic_variables, expected);
// }
//
// TEST(DegenerateSimplexMethodTests, ZerosInb2) {
//   size_t N = 10;
//
//   Matrix<Rational> A(2 * N, 3 * N, 0);
//
//   for (size_t i = 0; i < N; ++i) {
//     A[i, i] = Rational(i + 1);
//     A[N + i, N + i] = -Rational(i + 1);
//   }
//
//   Matrix<Rational> b(2 * N, 1, 0);
//   Matrix<Rational> c(1, 3 * N, 1);
//
//   auto bfs = BoundedSimplexMethod(CSCMatrix(A), b, c).find_bfs();
//
//   ASSERT_TRUE(bfs.has_value());
//
//   std::vector<size_t> expected(2 * N);
//   std::iota(expected.begin(), expected.end(), 0);
//
//   ASSERT_SETS_EQ(bfs->basic_variables, expected);
// }
