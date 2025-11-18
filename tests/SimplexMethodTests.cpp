#include <gtest/gtest.h>

#include <iostream>
#include <variant>

#include "linear/BigInteger.h"
#include "linear/Matrix.h"
#include "linear/SimplexMethod.h"

TEST(SimplexMethodTests, TableauInitialization) {
  {
    Matrix<Rational> A = {
        {3, 1, 1, 0, 0},
        {1, -1, 0, 1, 0},
        {0, 1, 0, 0, 1},
    };

    Matrix<Rational> b = {{6}, {2}, {3}};
    Matrix<Rational> c = {{-2, -1, 0, 0, 0}};

    Matrix<Rational> bfs = {{0}, {0}, {6}, {2}, {3}};

    SimplexMethod solver(A, b, c);
    auto tableau =
        solver.initialize_tableau(BFS<Rational>::construct_nondegenerate(bfs));

    Matrix<Rational> expected = {{6, 3, 1, 1, 0, 0},
                                 {2, 1, -1, 0, 1, 0},
                                 {3, 0, 1, 0, 0, 1},
                                 {0, 2, 1, 0, 0, 0}};

    ASSERT_EQ(expected, tableau);
  }

  {
    Matrix<Rational> A = {
        {1, -1, 1, 0},
        {2, 1, 0, 1},
    };

    Matrix<Rational> b = {{1}, {3}};
    Matrix<Rational> c = {{2, 1, 1, -1}};

    Matrix<Rational> bfs = {{0}, {0}, {1}, {3}};

    SimplexMethod solver(A, b, c);
    auto tableau =
        solver.initialize_tableau(BFS<Rational>::construct_nondegenerate(bfs));

    Matrix<Rational> expected = {
        {1, 1, -1, 1, 0}, {3, 2, 1, 0, 1}, {-2, -3, -3, 0, 0}};

    ASSERT_EQ(expected, tableau);
  }
}

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

TEST(SimplexMethodTests, FullSimple) {
  {
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

  {
    Matrix<Rational> A = {
        {1, 1, -1, 1},
        {1, 14, 10, -10},
    };

    Matrix<Rational> b = {{2}, {24}};
    Matrix<Rational> c = {{1, 2, 3, -4}};

    Matrix<Rational> bfs = {
        {Rational{0}}, {Rational{11} / 6}, {Rational{0}}, {Rational{1} / 6}};

    SimplexMethod solver(A, b, c);
    auto solution =
        solver.solve_from(BFS<Rational>::construct_nondegenerate(bfs));

    Matrix<Rational> expected = {{4}, {0}, {2}, {0}};

    ASSERT_EQ(std::get<FiniteLPSolution<Rational>>(solution).point, expected);
    ASSERT_EQ(std::get<FiniteLPSolution<Rational>>(solution).value, 10);
  }

  {
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

  {
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
