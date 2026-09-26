#include <gtest/gtest.h>
#include <gmpxx.h>

#include "support/GMPRational.h"

using Traits = FieldTraits<GMPRational>;

TEST(GMPTest, ItWorks) {
  mpz_class x = 123;

  x += 321;

  ASSERT_EQ(x, 444);
}

TEST(GMPTest, Arithmetic) {
  const GMPRational third = GMPRational(1) / 3;

  ASSERT_EQ(third + third + third, 1);
  ASSERT_EQ(-third * 3, -1);
  ASSERT_LT(third, GMPRational(1) / 2);
  ASSERT_EQ(abs(-third), third);
  ASSERT_EQ(to_string(-third), "-1/3");
  ASSERT_EQ(std::format("{}", third), "1/3");
}

TEST(GMPTest, FloorAndFractional) {
  ASSERT_EQ(Traits::floor(GMPRational(7) / 2), 3);
  ASSERT_EQ(Traits::floor(GMPRational(-7) / 2), -4);
  ASSERT_EQ(Traits::floor(5), 5);
  ASSERT_EQ(Traits::fractional(GMPRational(-7) / 2), GMPRational(1) / 2);
}

TEST(GMPTest, Exp2) {
  ASSERT_EQ(Traits::exp2(0), 1);
  ASSERT_EQ(Traits::exp2(10), 1024);
  ASSERT_EQ(Traits::exp2(-3), GMPRational(1) / 8);
}

TEST(GMPTest, FromString) {
  ASSERT_EQ(Traits::from_string("42"), 42);
  ASSERT_EQ(Traits::from_string("+42"), 42);
  ASSERT_EQ(Traits::from_string("-1.25"), GMPRational(-5) / 4);
  ASSERT_EQ(Traits::from_string(".5"), GMPRational(1) / 2);
  ASSERT_EQ(Traits::from_string("3."), 3);
  ASSERT_EQ(Traits::from_string("1.5e3"), 1500);
  ASSERT_EQ(Traits::from_string("25E-2"), GMPRational(1) / 4);
  ASSERT_EQ(Traits::from_string("0.1"), GMPRational(1) / 10);

  ASSERT_EQ(Traits::from_string(""), std::nullopt);
  ASSERT_EQ(Traits::from_string("."), std::nullopt);
  ASSERT_EQ(Traits::from_string("1.2.3"), std::nullopt);
  ASSERT_EQ(Traits::from_string("1e"), std::nullopt);
  ASSERT_EQ(Traits::from_string("1e+-2"), std::nullopt);
  ASSERT_EQ(Traits::from_string("abc"), std::nullopt);
}
