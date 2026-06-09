#pragma once

#include <ranges>

#include "linalg/Concepts.h"

class MatrixSpy {
 public:
  using FieldType = int;

  static size_t copy_constructor_calls;
  static size_t move_constructor_calls;
  static size_t copy_assignment_calls;
  static size_t move_assignment_calls;

  MatrixSpy() = default;

  MatrixSpy(const MatrixSpy&) { ++copy_constructor_calls; }
  MatrixSpy(MatrixSpy&&) noexcept { ++move_constructor_calls; }

  MatrixSpy& operator=(const MatrixSpy&) {
    ++copy_assignment_calls;
    return *this;
  }
  MatrixSpy& operator=(MatrixSpy&&) noexcept {
    ++move_assignment_calls;
    return *this;
  }

  static void reset_counters() {
    copy_constructor_calls = 0;
    move_constructor_calls = 0;
    copy_assignment_calls = 0;
    move_assignment_calls = 0;
  }

  auto entries() const {
    return std::ranges::empty_view<std::tuple<size_t, size_t, int>>{};
  }

  auto row_entries(size_t row) const {
    return std::ranges::empty_view<std::pair<size_t, int>>{};
  }

  auto col_entries(size_t col) const {
    return std::ranges::empty_view<std::pair<size_t, int>>{};
  }

  int operator[](size_t i, size_t j) const { return 0; }

  size_t cols() const { return 123; }
  size_t rows() const { return 123; }
  std::pair<size_t, size_t> shape() const { return {rows(), cols()}; }
};

inline size_t MatrixSpy::copy_constructor_calls = 0;
inline size_t MatrixSpy::move_constructor_calls = 0;
inline size_t MatrixSpy::copy_assignment_calls = 0;
inline size_t MatrixSpy::move_assignment_calls = 0;

static_assert(linalg::MatrixRange<MatrixSpy>);
static_assert(linalg::ColWiseMatrixRange<MatrixSpy>);
static_assert(linalg::RowWiseMatrixRange<MatrixSpy>);
static_assert(linalg::ElementWiseMatrixRange<MatrixSpy>);
