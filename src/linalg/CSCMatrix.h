#pragma once

#include <algorithm>
#include <ranges>
#include <vector>
#include <format>

#include "expr/SubColsExpr.h"
#include "field/FieldTraits.h"

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

  template <ElementWiseMatrixRange R>
    requires std::same_as<MatrixFieldType<R>, Field>
  explicit CSCMatrix(R&& matrix,
                     Field drop_tolerance = FieldTraits<Field>::tolerance)
      : CSCMatrix(matrix.rows(), 0) {
    using std::abs;

    for (size_t col = 0; col < matrix.cols(); ++col) {
      add_column();

      for (size_t row = 0; row < matrix.rows(); ++row) {
        if (abs(matrix[row, col]) > drop_tolerance) {
          push_to_last_column(row, matrix[row, col]);
        }
      }
    }
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
  template <typename F>
  void col_entries(size_t col, F&& f) const {
    for (size_t i = index_pointers_[col]; i < index_pointers_[col + 1]; ++i) {
      const auto [row, value] = entries_[i];
      f(row, value);
    }
  }

  template <typename F>
  void entries(F&& f) const {
    for (size_t col = 0; col < cols(); ++col) {
      for (size_t i = index_pointers_[col]; i < index_pointers_[col + 1]; ++i) {
        const auto [row, value] = entries_[i];

        f(row, col, value);
      }
    }
  }

  std::span<std::pair<size_t, Field>> get_column(size_t col) {
    return std::span{entries_.begin() + index_pointers_[col],
                     entries_.begin() + index_pointers_[col + 1]};
  }

  std::span<const std::pair<size_t, Field>> get_column(size_t col) const {
    return std::span{entries_.begin() + index_pointers_[col],
                     entries_.begin() + index_pointers_[col + 1]};
  }

  auto get_column_as_matrix(size_t col) const {
    return detail::SubColsExpr(*this, std::ranges::single_view{col});
  }

  template <IndicesRange R>
  auto select_columns(R&& cols) const {
    return detail::SubColsExpr(*this, std::forward<R>(cols));
  }

  // Has O(n) complexity, where n is the entries count in the column.
  std::optional<Field> at(size_t row, size_t col) const {
    Field sum = 0;
    bool found = false;

    for (const auto [other_row, value] : get_column(col)) {
      if (row == other_row) {
        sum += value;
        found = true;
      }
    }

    return found ? std::optional{sum} : std::nullopt;
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

  // Removes all entries and columns from the matrix. If matrix shape was
  // (rows, cols), then after clear it would be (rows, 0).
  void clear() {
    entries_.clear();
    index_pointers_.clear();

    index_pointers_.push_back(0);
  }
};

template <ElementWiseMatrixRange R>
CSCMatrix(R&& matrix, MatrixFieldType<R> drop_tolerance =
                          FieldTraits<MatrixFieldType<R>>::tolerance)
    -> CSCMatrix<MatrixFieldType<R>>;

}  // namespace linalg
