#pragma once

#include <set>

void assert_sets_equal(auto&& left_range, auto&& right_range) {
  auto left_set =
      std::set(std::from_range, std::forward<decltype(left_range)>(left_range));

  auto right_set = std::set(std::from_range,
                            std::forward<decltype(right_range)>(right_range));

  ASSERT_EQ(left_set, right_set);
}

#define ASSERT_SETS_EQ(left, right) assert_sets_equal(left, right)
