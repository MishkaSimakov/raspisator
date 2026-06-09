#pragma once

#include <algorithm>
#include <ranges>
#include <vector>

#include "expr/SubColsExpr.h"

#include "Arithmetics.h"

namespace linalg {

// Stores sparse matrix in Compacted Sparse Column form. For each column it
// stores a list of entries (row + value) as a contiguous subrange in entries_.
// In one column entries are:
// 1. Unordered
// 2. Possibly duplicated (same row may appear twice). In this case values are
// added up.
// 3. Bounded by matrix height (row >= rows() is invalid).
template <typename Field>
class CSCMatrix {
  std::vector<std::pair<size_t, Field>> entries_;
  std::vector<size_t> index_pointers_;

  size_t rows_cnt_;

  explicit CSCMatrix(size_t rows, size_t cols)
      : index_pointers_(cols + 1, 0), rows_cnt_(rows) {}

 public:
  using FieldType = Field;

  CSCMatrix() : CSCMatrix(0, 0) {}

  static CSCMatrix zeros(size_t rows = 0, size_t cols = 0) {
    return CSCMatrix(rows, cols);
  }

  template <typename R>
    requires std::ranges::range<R> &&
             std::convertible_to<std::ranges::range_value_t<R>,
                                 std::pair<size_t, Field>>
  void add_column(R&& column) {
    entries_.insert(entries_.end(), std::ranges::begin(column),
                    std::ranges::end(column));
    index_pointers_.push_back(entries_.size());

    const size_t begin = index_pointers_[index_pointers_.size() - 2];
    const size_t end = index_pointers_[index_pointers_.size() - 1];

    for (size_t i = begin; i < end; ++i) {
      if (entries_[i].first >= rows()) {
        throw std::invalid_argument(
            std::format("Row index {} is invalid for matrix with height {}.",
                        entries_[i].first, rows()));
      }
    }
  }

  void add_column() { index_pointers_.push_back(index_pointers_.back()); }

  void push_to_last_column(size_t row, Field value) {
    if (row >= rows()) {
      throw std::out_of_range(std::format(
          "Row index {} is invalid for matrix with height {}.", row, rows()));
    }
    if (cols() == 0) {
      throw std::out_of_range(
          "Can't push to last column because there are no columns.");
    }

    entries_.emplace_back(row, value);
    ++index_pointers_.back();
  }

  //
  auto entries() const {
    return std::views::iota(size_t{0}, index_pointers_.size() - 1) |
           std::views::transform([this](size_t col) {
             auto begin = index_pointers_[col];
             auto end = index_pointers_[col + 1];

             return std::views::iota(begin, end) |
                    std::views::transform([this, col](size_t idx) {
                      const auto [row, value] = entries_[idx];
                      return std::tuple{row, col, value};
                    });
           }) |
           std::views::join;
  }

  template <IndicesRange R>
  auto get_columns(R&& cols) const {
    return detail::SubColsExpr(*this, std::forward<R>(cols));
  }

  auto get_column(size_t col) const {
    return get_columns(std::views::single(col));
  }

  std::span<std::pair<size_t, Field>> col_entries(size_t col) {
    return {entries_.begin() + index_pointers_[col],
            entries_.begin() + index_pointers_[col + 1]};
  }

  std::span<const std::pair<size_t, Field>> col_entries(size_t col) const {
    return {entries_.begin() + index_pointers_[col],
            entries_.begin() + index_pointers_[col + 1]};
  }

  //
  size_t rows() const { return rows_cnt_; }
  size_t cols() const { return index_pointers_.size() - 1; }
  std::pair<size_t, size_t> shape() const { return {rows(), cols()}; }

  size_t entries_count() const { return entries_.size(); }

  //
  void resize(size_t new_rows, size_t new_cols) {
    index_pointers_.resize(new_cols + 1, index_pointers_.back());

    // resize rows count
    if (new_rows >= rows()) {
      entries_.resize(index_pointers_.back());
      rows_cnt_ = new_rows;
      return;
    }

    size_t offset = 0;
    for (size_t col = 0; col < new_cols; ++col) {
      const size_t column_start = index_pointers_[col] + offset;

      for (size_t i = column_start; i < index_pointers_[col + 1]; ++i) {
        if (entries_[i].first < new_rows) {
          entries_[i - offset] = entries_[i];
        } else {
          ++offset;
        }
      }

      index_pointers_[col + 1] -= offset;
    }

    entries_.resize(index_pointers_.back());
    rows_cnt_ = new_rows;
  }

  void map_rows(std::span<const size_t> map) {
    // TODO: dimensions check
    for (auto& [row, value] : entries_) {
      row = map[row];
    }
  }
};

}  // namespace linalg
