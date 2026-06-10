#pragma once

#include <map>
#include <tuple>

#include "linalg/Concepts.h"

template <linalg::MatrixRange L, typename R>
  requires std::ranges::range<R> &&
           std::same_as<std::ranges::range_value_t<R>,
                        std::tuple<size_t, size_t, linalg::MatrixFieldType<L>>>
bool check_entries(L&& left, R&& right) {
  using Field = linalg::MatrixFieldType<L>;

  std::map<std::pair<size_t, size_t>, Field> left_map;
  left.entries([&](size_t i, size_t j, Field value) {
    left_map[std::pair{i, j}] += value;
  });

  std::map<std::pair<size_t, size_t>, Field> right_map;
  for (const auto [i, j, value] : right) {
    right_map[std::pair{i, j}] += value;
  }

  return left_map == right_map;
}

template <linalg::ColWiseMatrixRange L, typename R>
  requires std::ranges::range<R> &&
           std::same_as<std::ranges::range_value_t<R>,
                        std::pair<size_t, linalg::MatrixFieldType<L>>>
bool check_col_entries(L&& left, size_t col, R&& right) {
  using Field = linalg::MatrixFieldType<L>;

  std::map<size_t, Field> left_map;
  left.col_entries(col, [&](size_t i, Field value) { left_map[i] += value; });

  std::map<size_t, Field> right_map;
  for (const auto [i, value] : right) {
    right_map[i] += value;
  }

  return left_map == right_map;
}

template <linalg::RowWiseMatrixRange L, typename R>
  requires std::ranges::range<R> &&
           std::same_as<std::ranges::range_value_t<R>,
                        std::pair<size_t, linalg::MatrixFieldType<L>>>
bool check_row_entries(L&& left, size_t row, R&& right) {
  using Field = linalg::MatrixFieldType<L>;

  std::map<size_t, Field> left_map;
  left.row_entries(row, [&](size_t i, Field value) { left_map[i] += value; });

  std::map<size_t, Field> right_map;
  for (const auto [i, value] : right) {
    right_map[i] += value;
  }

  return left_map == right_map;
}

#define ASSERT_ENTRIES_EQ(matrix, expected_entries) \
  ASSERT_TRUE(check_entries(matrix, expected_entries))

#define ASSERT_COL_ENTRIES_EQ(matrix, col, expected_entries) \
  ASSERT_TRUE(check_col_entries(matrix, col, expected_entries))

#define ASSERT_ROW_ENTRIES_EQ(matrix, row, expected_entries) \
  ASSERT_TRUE(check_row_entries(matrix, row, expected_entries))
