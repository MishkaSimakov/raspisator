#include <gtest/gtest.h>

#include <sstream>

#include "linalg/Linalg.h"

TEST(PrintTest, MatrixRange) {
  Matrix<int> A = {{1}};
  Matrix<int> B = {{2}};

  std::stringstream ss;
  ss << (A + B);

  ASSERT_TRUE(ss.str().contains("3"));
}

TEST(PrintTest, Dense) {
  Matrix<int> A = {{1, 2}, {3, 4}};

  std::stringstream ss;
  ss << A;

  ASSERT_TRUE(ss.str().contains("1"));
  ASSERT_TRUE(ss.str().contains("2"));
  ASSERT_TRUE(ss.str().contains("3"));
  ASSERT_TRUE(ss.str().contains("4"));
}
