#include <gtest/gtest.h>

#include <chrono>

#include "linear/BigInteger.h"
#include "utils/Accumulators.h"

TEST(ArithmeticMeanTests, TestEmpty) {
  ArithmeticMean<Rational> mean;

  ASSERT_EQ(mean.count(), 0);
  ASSERT_EQ(mean.sum(), 0);
  ASSERT_EQ(mean.get(), std::nullopt);
}

TEST(ArithmeticMeanTests, TestCorrectCasting1) {
  ArithmeticMean<Rational> mean;

  mean.record(1);
  mean.record(2);

  ASSERT_EQ(*mean, Rational{3} / 2);
}

TEST(ArithmeticMeanTests, TestCorrectCasting2) {
  ArithmeticMean<double> mean;

  mean.record(1);
  mean.record(2);

  ASSERT_DOUBLE_EQ(*mean, 1.5);
}

TEST(ArithmeticMeanTests, TestCorrectCasting3) {
  using namespace std::chrono_literals;

  ArithmeticMean<std::chrono::seconds> mean;

  mean.record(5s);
  mean.record(15s);

  ASSERT_EQ(*mean, 10s);
}
