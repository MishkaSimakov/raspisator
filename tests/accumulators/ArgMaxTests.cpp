#include <gtest/gtest.h>

#include "utils/Accumulators.h"

TEST(ArgMaxTests, EmptyArgMax) {
  ArgMaximum<int> max;

  ASSERT_EQ(max.get(), std::nullopt);
  ASSERT_EQ(max.has_value(), false);
}

TEST(ArgMaxTests, SimpleTest) {
  ArgMaximum<int> max;

  max.record(1, 2);
  max.record(2, -2);
  max.record(3, 4);

  ASSERT_EQ(max->index, 3);
  ASSERT_EQ(max->max, 4);
}

TEST(ArgMaxTests, SelectsSmallestIndex) {
  ArgMaximum<int> max;

  max.record(1, 2);
  max.record(4, 6);
  max.record(2, 6);
  max.record(3, 6);
  max.record(5, -10);

  ASSERT_EQ(max->index, 2);
  ASSERT_EQ(max->max, 6);
}
