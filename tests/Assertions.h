#pragma once

#include <set>

void assert_sets_equal(auto&& left_range, auto&& right_range) {
  std::set left_set(left_range.begin(), left_range.end());
  std::set right_set(right_range.begin(), right_range.end());

  ASSERT_EQ(left_set, right_set);
}

#define ASSERT_SETS_EQ(left, right) assert_sets_equal(left, right)
