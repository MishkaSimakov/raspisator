#pragma once

#include <map>

#include "linalg/Concepts.h"

// TODO: Field default constructor must initialize it with 0
template <typename L, typename R>
  requires linalg::DoublesRange<L, linalg::DoublesRangeFieldType<L>> &&
           linalg::DoublesRange<R, linalg::DoublesRangeFieldType<R>> &&
           std::same_as<linalg::DoublesRangeFieldType<L>,
                        linalg::DoublesRangeFieldType<R>>
bool check_doubles_ranges_eq(L&& left, R&& right) {
  using Field = linalg::DoublesRangeFieldType<L>;

  std::map<size_t, Field> left_map;
  for (const auto& [i, value] : left) {
    left_map[i] += value;
  }

  std::map<size_t, Field> right_map;
  for (const auto& [i, value] : right) {
    right_map[i] += value;
  }

  return left_map == right_map;
}

template <typename L, typename R>
  requires linalg::TriplesRange<L, linalg::TriplesRangeFieldType<L>> &&
           linalg::TriplesRange<R, linalg::TriplesRangeFieldType<R>> &&
           std::same_as<linalg::TriplesRangeFieldType<L>,
                        linalg::TriplesRangeFieldType<R>>
bool check_triples_ranges_eq(L&& left, R&& right) {
  using Field = linalg::TriplesRangeFieldType<L>;

  std::map<std::pair<size_t, size_t>, Field> left_map;
  for (const auto& [i, j, value] : left) {
    left_map[std::pair{i, j}] += value;
  }

  std::map<std::pair<size_t, size_t>, Field> right_map;
  for (const auto& [i, j, value] : right) {
    right_map[std::pair{i, j}] += value;
  }

  return left_map == right_map;
}

#define ASSERT_DOUBLES_RANGES_EQ(left, right) \
  ASSERT_TRUE(check_doubles_ranges_eq(left, right))

#define ASSERT_TRIPLES_RANGES_EQ(left, right) \
  ASSERT_TRUE(check_triples_ranges_eq(left, right))
