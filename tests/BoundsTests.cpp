#include <gtest/gtest.h>

#include "linear/BigInteger.h"
#include "linear/model/Bound.h"

TEST(BoundsTests, SimpleMultiplication) {
  Bound<Rational> bound(-10, 10);

  bound *= 5;

  ASSERT_TRUE(bound.lower == -50 && bound.upper == 50);
}

TEST(BoundsTests, NegativeMultiplication) {
  Bound<Rational> bound(-5, 10);

  bound *= -5;

  ASSERT_TRUE(bound.lower == -50 && bound.upper == 25);
}

TEST(BoundsTests, UnboundedMultiplication) {
  Bound<Rational> bound(-5, std::nullopt);

  bound *= 5;

  ASSERT_TRUE(bound.lower == -25 && bound.upper == std::nullopt);
}

TEST(BoundsTests, Subtraction) {
  Bound<Rational> left(-5, 10);
  Bound<Rational> right(10, 20);

  auto diff = left - right;

  ASSERT_TRUE(diff.lower == -25 && diff.upper == 0);
}

TEST(BoundsTests, Intersection1) {
  Bound<Rational> left(-5, 10);
  Bound<Rational> right(5, 20);

  auto intersection = left ^ right;

  ASSERT_TRUE(intersection.lower == 5 && intersection.upper == 10);
}

TEST(BoundsTests, Intersection2) {
  Bound<Rational> left(-10, 10);
  Bound<Rational> right(-10, 5);

  auto intersection = left ^ right;

  ASSERT_TRUE(intersection.lower == -10 && intersection.upper == 5);
}

TEST(BoundsTests, Intersection3) {
  Bound<Rational> left(-10, 10);
  Bound<Rational> right(-10, 10);

  auto intersection = left ^ right;

  ASSERT_TRUE(intersection.lower == -10 && intersection.upper == 10);
}

TEST(BoundsTests, Intersection4) {
  Bound<Rational> left(-10, 10);
  Bound<Rational> right(std::nullopt, std::nullopt);

  auto intersection = left ^ right;

  ASSERT_TRUE(intersection.lower == -10 && intersection.upper == 10);
}

TEST(BoundsTests, Intersection5) {
  Bound<Rational> left(0, 1);
  Bound<Rational> right(-5, 6);

  auto intersection = left ^ right;

  ASSERT_TRUE(intersection.lower == 0 && intersection.upper == 1);
}

TEST(BoundsTests, Intersection6) {
  Bound<Rational> left(0, 1);
  Bound<Rational> right(2, 3);

  auto intersection = left ^ right;

  ASSERT_TRUE(intersection.is_infeasible());
}
