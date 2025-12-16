#include <gtest/gtest.h>

#include "utils/Accumulators.h"

TEST(ArgMinTests, EmptyArgMin) {
  ArgMinimum<int> min;

  ASSERT_EQ(min.argmin(), std::nullopt);
  ASSERT_EQ(min.min(), std::nullopt);
}

TEST(ArgMinTests, SimpleTest) {
  ArgMinimum<int> min;

  min.record(1, 2);
  min.record(2, -2);
  min.record(3, 3);

  ASSERT_EQ(min.argmin(), 2);
  ASSERT_EQ(min.min(), -2);
}

TEST(ArgMinTests, SelectsSmallestIndex) {
  ArgMinimum<int> min;

  min.record(1, 2);
  min.record(4, -2);
  min.record(2, -2);
  min.record(3, -2);
  min.record(5, 10);

  ASSERT_EQ(min.argmin(), 2);
  ASSERT_EQ(min.min(), -2);
}
