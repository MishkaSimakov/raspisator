#include <gtest/gtest.h>

#include "problem/Bound.h"
#include "support/GMPRational.h"

TEST(BoundsTests, SimpleMultiplication) {
  Bound<GMPRational> bound(-10, 10);

  bound *= 5;

  ASSERT_TRUE(bound.lower == -50 && bound.upper == 50);
}

TEST(BoundsTests, NegativeMultiplication) {
  Bound<GMPRational> bound(-5, 10);

  bound *= -5;

  ASSERT_TRUE(bound.lower == -50 && bound.upper == 25);
}

TEST(BoundsTests, UnboundedMultiplication) {
  Bound<GMPRational> bound(-5, std::nullopt);

  bound *= 5;

  ASSERT_TRUE(bound.lower == -25 && bound.upper == std::nullopt);
}

TEST(BoundsTests, Subtraction) {
  Bound<GMPRational> left(-5, 10);
  Bound<GMPRational> right(10, 20);

  auto diff = left - right;

  ASSERT_TRUE(diff.lower == -25 && diff.upper == 0);
}

TEST(BoundsTests, Intersection1) {
  Bound<GMPRational> left(-5, 10);
  Bound<GMPRational> right(5, 20);

  auto intersection = left ^ right;

  ASSERT_TRUE(intersection.lower == 5 && intersection.upper == 10);
}

TEST(BoundsTests, Intersection2) {
  Bound<GMPRational> left(-10, 10);
  Bound<GMPRational> right(-10, 5);

  auto intersection = left ^ right;

  ASSERT_TRUE(intersection.lower == -10 && intersection.upper == 5);
}

TEST(BoundsTests, Intersection3) {
  Bound<GMPRational> left(-10, 10);
  Bound<GMPRational> right(-10, 10);

  auto intersection = left ^ right;

  ASSERT_TRUE(intersection.lower == -10 && intersection.upper == 10);
}

TEST(BoundsTests, Intersection4) {
  Bound<GMPRational> left(-10, 10);
  Bound<GMPRational> right(std::nullopt, std::nullopt);

  auto intersection = left ^ right;

  ASSERT_TRUE(intersection.lower == -10 && intersection.upper == 10);
}

TEST(BoundsTests, Intersection5) {
  Bound<GMPRational> left(0, 1);
  Bound<GMPRational> right(-5, 6);

  auto intersection = left ^ right;

  ASSERT_TRUE(intersection.lower == 0 && intersection.upper == 1);
}

TEST(BoundsTests, Intersection6) {
  Bound<GMPRational> left(0, 1);
  Bound<GMPRational> right(2, 3);

  auto intersection = left ^ right;

  ASSERT_TRUE(intersection.is_infeasible());
}
