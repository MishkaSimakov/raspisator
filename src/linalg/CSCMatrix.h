#pragma once

#include <algorithm>
#include <ranges>
#include <vector>

#include "expr/SubColsExpr.h"
#include "expr/TransposedExpr.h"

#include "Arithmetics.h"

namespace linalg {

template <typename Field>
class CSCMatrix {
  std::vector<std::pair<size_t, Field>> entries_;
  std::vector<size_t> index_pointers_;

  size_t rows_cnt_;

  explicit CSCMatrix(size_t rows, size_t cols)
      : index_pointers_(cols + 1, 0), rows_cnt_(rows) {}

  // performs binary search on entries_
  std::optional<size_t> get_entry_index(size_t row, size_t col) const {
    size_t left = index_pointers_[col];
    size_t right = index_pointers_[col + 1];

    while (left + 1 < right) {
      const size_t middle = (left + right) / 2;

      if (entries_[middle].first == row) {
        return middle;
      }

      if (entries_[middle].first < row) {
        left = middle;
      } else {
        right = middle;
      }
    }

    return entries_[left].first == row ? std::optional{left} : std::nullopt;
  }

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

    if (begin == end) {
      return;
    }

    for (size_t i = begin; i < end; ++i) {
      if (entries_[i].first >= rows()) {
        throw std::invalid_argument(
            std::format("Row index {} is invalid for matrix with height {}.",
                        entries_[i].first, rows()));
      }
    }

    std::ranges::sort(
        entries_.begin() + begin, entries_.begin() + end, {},
        [](const std::pair<size_t, Field>& entry) { return entry.first; });

    size_t shift = 0;
    size_t row = entries_[begin].first;

    for (size_t i = begin + 1; i < end; ++i) {
      if (entries_[i].first == row) {
        ++shift;
        entries_[i - shift].second += entries_[i].second;
      } else {
        entries_[i - shift] = entries_[i];
        row = entries_[i].first;
      }
    }

    entries_.resize(entries_.size() - shift);
    index_pointers_.back() -= shift;
  }

  void add_column() { index_pointers_.push_back(index_pointers_.back()); }

  //
  Field operator[](size_t row, size_t col) const {
    const auto index = get_entry_index(row, col);

    return index ? entries_[*index].second : 0;
  }

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

  auto transposed() { return detail::TransposedExpr(*this); }
  auto transposed() const { return detail::TransposedExpr(*this); }

  void map_rows(std::span<const size_t> map) {
    // TODO: dimensions check
    for (auto& [row, value] : entries_) {
      row = map[row];
    }
  }
};

}  // namespace linalg
