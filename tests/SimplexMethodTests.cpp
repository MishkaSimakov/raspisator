#include <gtest/gtest.h>

#include <iostream>
#include <variant>

#include "linear/BigInteger.h"
#include "linear/CheckBFS.h"
#include "linear/SimplexMethod.h"
#include "linear/matrix/Matrix.h"

TEST(SimplexMethodTests, SimplexMethodStartingInSolution) {
  Matrix<Rational> A = {
      {1, -1, 1, 0},
      {2, 1, 0, 1},
  };

  Matrix<Rational> b = {{1}, {3}};
  Matrix<Rational> c = {{2, 1, 1, -1}};

  Matrix<Rational> bfs = {{0}, {3}, {4}, {0}};

  SimplexMethod solver(A, b, c);
  auto solution =
      solver.solve_from(BFS<Rational>::construct_nondegenerate(bfs));

  ASSERT_EQ(std::get<FiniteLPSolution<Rational>>(solution).point, bfs);
  ASSERT_EQ(std::get<FiniteLPSolution<Rational>>(solution).value, 7);
}

TEST(SimplexMethodTests, FullSimple1) {
  Matrix<Rational> A = {
      {1, -1, 1, 0},
      {2, 1, 0, 1},
  };

  Matrix<Rational> b = {{1}, {3}};
  Matrix<Rational> c = {{2, 1, 1, -1}};

  Matrix<Rational> bfs = {{0}, {0}, {1}, {3}};

  SimplexMethod solver(A, b, c);
  auto solution =
      solver.solve_from(BFS<Rational>::construct_nondegenerate(bfs));

  Matrix<Rational> expected = {{0}, {3}, {4}, {0}};

  ASSERT_EQ(std::get<FiniteLPSolution<Rational>>(solution).point, expected);
  ASSERT_EQ(std::get<FiniteLPSolution<Rational>>(solution).value, 7);
}

TEST(SimplexMethodTests, FullSimple2) {
  Matrix<Rational> A = {
      {1, 1, -1, 1},
      {1, 14, 10, -10},
  };

  Matrix<Rational> b = {{2}, {24}};
  Matrix<Rational> c = {{1, 2, 3, -4}};

  Matrix bfs = {
      {Rational{0}}, {Rational{11} / 6}, {Rational{0}}, {Rational{1} / 6}};

  SimplexMethod solver(A, b, c);
  auto solution =
      solver.solve_from(BFS<Rational>::construct_nondegenerate(bfs));

  Matrix<Rational> expected = {{4}, {0}, {2}, {0}};

  ASSERT_EQ(std::get<FiniteLPSolution<Rational>>(solution).point, expected);
  ASSERT_EQ(std::get<FiniteLPSolution<Rational>>(solution).value, 10);
}

TEST(SimplexMethodTests, FullSimple3) {
  Matrix<Rational> A = {{1, 1}};
  Matrix<Rational> b = {{1}};
  Matrix<Rational> c = {{1, 2}};

  Matrix<Rational> bfs = {{1}, {0}};

  SimplexMethod solver(A, b, c);
  auto solution =
      solver.solve_from(BFS<Rational>::construct_nondegenerate(bfs));

  Matrix<Rational> expected = {{0}, {1}};

  ASSERT_EQ(std::get<FiniteLPSolution<Rational>>(solution).point, expected);
  ASSERT_EQ(std::get<FiniteLPSolution<Rational>>(solution).value, 2);
}

TEST(SimplexMethodTests, FullSimple4) {
  Matrix<Rational> A = {{1, 1}};
  Matrix<Rational> b = {{1}};
  Matrix<Rational> c = {{2, 1}};

  Matrix<Rational> bfs = {{0}, {1}};

  SimplexMethod solver(A, b, c);
  auto solution =
      solver.solve_from(BFS<Rational>::construct_nondegenerate(bfs));

  Matrix<Rational> expected = {{1}, {0}};

  ASSERT_EQ(std::get<FiniteLPSolution<Rational>>(solution).point, expected);
  ASSERT_EQ(std::get<FiniteLPSolution<Rational>>(solution).value, 2);
}

TEST(SimplexMethodTests, InfiniteSolutionDetection) {
  {
    Matrix<Rational> A = {{1, 0}};
    Matrix<Rational> b = {{1}};
    Matrix<Rational> c = {{0, 1}};

    Matrix<Rational> bfs = {{1}, {0}};

    SimplexMethod solver(A, b, c);
    auto solution =
        solver.solve_from(BFS<Rational>::construct_nondegenerate(bfs));

    ASSERT_TRUE(std::holds_alternative<InfiniteSolution>(solution));
  }
}

TEST(SimplexMethodTests, FullWithFindingBFS) {
  Matrix<Rational> A = {
      {1, 1, -1, 1},
      {1, 14, 10, -10},
  };

  Matrix<Rational> b = {{2}, {24}};
  Matrix<Rational> c = {{1, 2, 3, -4}};

  SimplexMethod solver(A, b, c);
  auto solution = solver.solve();

  Matrix<Rational> expected = {{4}, {0}, {2}, {0}};

  ASSERT_EQ(std::get<FiniteLPSolution<Rational>>(solution).point, expected);
  ASSERT_EQ(std::get<FiniteLPSolution<Rational>>(solution).value, 10);
}

TEST(SimplexMethodTests, ReconstructBFS) {
  Matrix<Rational> A = {
      {1, 1, -1, 1, 0, 0},
      {1, 14, 10, -15, 0, 0},
      {0, 0, 0, 0, 1, -1},
  };

  Matrix<Rational> c = {{1, 2, 3, -4, 0, 0}};

  // all combinations of variables
  for (size_t i1 = 0; i1 <= 1; ++i1) {
    for (size_t i2 = 0; i2 <= 1; ++i2) {
      for (size_t i3 = 0; i3 <= 1; ++i3) {
        for (size_t i4 = 0; i4 <= 1; ++i4) {
          if (i1 + i2 + i3 + i4 != 2) {
            continue;
          }

          auto old_bfs = BFS<Rational>::construct_nondegenerate(
              Matrix<Rational>{{i1}, {i2}, {i3}, {i4}, {-10}, {0}});

          auto b = A * old_bfs.point;

          auto bfs = SimplexMethod(A, b, c).reconstruct_bfs(old_bfs, 4);

          ASSERT_TRUE(bfs.has_value());
          ASSERT_TRUE(check_bfs(A, b, *bfs));
        }
      }
    }
  }
}
