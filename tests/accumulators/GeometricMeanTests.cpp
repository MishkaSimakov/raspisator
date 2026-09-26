#include <gtest/gtest.h>

#include <limits>

#include "utils/Accumulators.h"

TEST(GeometricMeanTests, TestEmpty) {
  GeometricMean<double> mean;

  ASSERT_EQ(mean.count(), 0);
  ASSERT_EQ(mean.has_value(), false);
  ASSERT_EQ(mean.get(), std::nullopt);
}

TEST(GeometricMeanTests, SimpleTest) {
  GeometricMean<double> mean;

  mean.record(2);
  mean.record(8);

  ASSERT_EQ(mean.count(), 2);
  ASSERT_DOUBLE_EQ(*mean, 4.0);
}

TEST(GeometricMeanTests, ThreeValues) {
  GeometricMean<double> mean;

  mean.record(1);
  mean.record(3);
  mean.record(9);

  ASSERT_NEAR(*mean, 3.0, 1e-12);
}

TEST(GeometricMeanTests, SingleValue) {
  GeometricMean<double> mean;

  mean.record(42);

  ASSERT_DOUBLE_EQ(*mean, 42.0);
}

TEST(GeometricMeanTests, AllOnes) {
  GeometricMean<double> mean;

  mean.record(1);
  mean.record(1);
  mean.record(1);

  ASSERT_DOUBLE_EQ(*mean, 1.0);
}

// The product of these values is far above std::numeric_limits<double>::max(),
// so any implementation accumulating the product directly would overflow to
// +inf. Accumulating in log-space keeps the result finite and correct.
TEST(GeometricMeanTests, NoOverflowOnLargeValues) {
  GeometricMean<double> mean;

  for (int i = 0; i < 100; ++i) {
    mean.record(1e300);
  }

  ASSERT_TRUE(std::isfinite(*mean));
  ASSERT_NEAR(*mean, 1e300, 1e300 * 1e-12);
}

// Symmetric to the overflow case: the product of many tiny values underflows
// to 0 in a naive implementation, while log-space accumulation preserves it.
TEST(GeometricMeanTests, NoUnderflowOnSmallValues) {
  GeometricMean<double> mean;

  for (int i = 0; i < 100; ++i) {
    mean.record(1e-300);
  }

  ASSERT_GT(*mean, 0.0);
  ASSERT_NEAR(*mean, 1e-300, 1e-300 * 1e-12);
}

// Mixing very large and very small factors: a naive product would overflow on
// the large factors before the small ones could bring it back into range.
TEST(GeometricMeanTests, NoOverflowOnMixedMagnitudes) {
  GeometricMean<double> mean;

  for (int i = 0; i < 50; ++i) {
    mean.record(1e300);
    mean.record(1e-300);
  }

  ASSERT_TRUE(std::isfinite(*mean));
  ASSERT_NEAR(*mean, 1.0, 1e-6);
}
