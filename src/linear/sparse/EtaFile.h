#pragma once

#include <type_traits>
#include <vector>

#include "CSCMatrix.h"
#include "linear/matrix/Matrix.h"

template <typename T, bool is_const>
using const_if = std::conditional_t<is_const, const T, T>;

namespace linalg {

template <typename Field>
class EtaFile {
  struct Entry {
    size_t begin;
    size_t index;
    bool is_column;
    bool is_removed;
  };

  template <bool is_const>
  struct EntryView {
    using ValueT = const_if<std::pair<size_t, Field>, is_const>;

    std::span<ValueT> values;

    const size_t index;
    const bool is_column;

    EntryView(std::span<ValueT> values, size_t index, bool is_column)
        : values(values), index(index), is_column(is_column) {}

    EntryView(EntryView<false> other)
        : values(other.values),
          index(other.index),
          is_column(other.is_column) {}
  };

  std::vector<std::pair<size_t, Field>> values_;
  std::vector<Entry> entries_;

  template <bool is_const>
  class iterator_base {
    using EntryT = const_if<Entry, is_const>;
    using ValueT = const_if<std::pair<size_t, Field>, is_const>;

    using ValuesT = const_if<std::vector<std::pair<size_t, Field>>, is_const>;

    EntryT* entry_;
    EntryT* end_;
    ValuesT* values_;

    iterator_base(EntryT* entry, EntryT* end, ValuesT& values)
        : entry_(entry), end_(end), values_(&values) {}

   public:
    using iterator_category = std::bidirectional_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = EntryView<is_const>;
    using reference = EntryView<is_const>;

    iterator_base() : entry_(nullptr), end_(nullptr), values_(nullptr) {}

    iterator_base(iterator_base<false> itr)
      requires(is_const)
        : entry_(itr.entry_), end_(itr.end_), values_(itr.values_) {}

    EntryView<is_const> operator*() const {
      size_t end = entry_ + 1 == end_ ? values_->size() : (entry_ + 1)->begin;

      return EntryView<is_const>(
          std::span{values_->begin() + entry_->begin, values_->begin() + end},
          entry_->index, entry_->is_column);
    }

    iterator_base& operator++() {
      ++entry_;
      while (entry_ != end_ && entry_->is_removed) {
        ++entry_;
      }

      return *this;
    }
    iterator_base operator++(int) {
      iterator_base tmp = *this;
      ++(*this);
      return tmp;
    }
    iterator_base& operator--() {
      --entry_;
      while (entry_->is_removed) {
        --entry_;
      }

      return *this;
    }
    iterator_base operator--(int) {
      iterator_base tmp = *this;
      --(*this);
      return tmp;
    }

    friend bool operator==(iterator_base left, iterator_base right) {
      return left.entry_ == right.entry_;
    }
    friend bool operator!=(iterator_base left, iterator_base right) {
      return left.entry_ != right.entry_;
    }

    friend EtaFile;
  };

 public:
  using iterator = iterator_base<false>;
  using const_iterator = iterator_base<true>;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

 private:
  static Matrix<Field> apply_impl(Matrix<Field> vector, EntryView<true> entry,
                                  bool transposed) {
    if (entry.is_column != transposed) {
      Field a = vector[entry.index, 0];
      vector[entry.index, 0] = 0;

      for (auto [row, value] : entry.values) {
        vector[row, 0] += a * value;
      }
    } else {
      Field dot = 0;

      for (auto [row, value] : entry.values) {
        dot += value * vector[row, 0];
      }

      vector[entry.index, 0] = dot;
    }

    return vector;
  }

 public:
  EtaFile() = default;

  void push_back(size_t column,
                 const std::vector<std::pair<size_t, Field>>& values) {
    entries_.push_back(Entry{
        .begin = values_.size(),
        .index = column,
        .is_column = true,
        .is_removed = false,
    });

    values_.append_range(values);
  }

  Matrix<Field> apply(Matrix<Field> vector, EntryView<true> entry) const {
    return apply_impl(std::move(vector), entry, false);
  }

  Matrix<Field> apply_transposed(Matrix<Field> vector,
                                 EntryView<true> entry) const {
    return apply_impl(std::move(vector), entry, true);
  }

  void purge() {
    // TODO: delete removed eta matrices and zeros
  }

  void clear() {
    values_.clear();
    entries_.clear();
  }

  size_t get_breadth() const {
    size_t result = 0;

    for (const auto& entry : entries_) {
      if (!entry.is_removed) {
        result = std::max(result, entry.index);
      }
    }

    for (auto [row, _] : values_) {
      result = std::max(result, row);
    }

    return result + 1;
  }

  iterator erase(iterator itr) {
    itr.entry_->is_removed = true;
    ++itr;

    return itr;
  }

  // iterators
  iterator begin() {
    Entry* end = entries_.data() + entries_.size();
    Entry* current = entries_.data();

    while (current != end && current->is_removed) {
      ++current;
    }

    return {current, end, values_};
  }
  const_iterator begin() const { return const_cast<EtaFile*>(this)->begin(); }
  const_iterator cbegin() const { return begin(); }

  iterator end() {
    Entry* end = entries_.data() + entries_.size();
    return {end, end, values_};
  }
  const_iterator end() const { return const_cast<EtaFile*>(this)->end(); }
  const_iterator cend() const { return end(); }

  reverse_iterator rbegin() { return std::make_reverse_iterator(end()); }
  const_reverse_iterator rbegin() const {
    return std::make_reverse_iterator(end());
  }
  const_reverse_iterator crbegin() const {
    return std::make_reverse_iterator(end());
  }

  reverse_iterator rend() { return std::make_reverse_iterator(begin()); }
  const_reverse_iterator rend() const {
    return std::make_reverse_iterator(begin());
  }
  const_reverse_iterator crend() const {
    return std::make_reverse_iterator(begin());
  }
};

}  // namespace linalg

template <typename Field>
std::ostream& operator<<(std::ostream& os, const linalg::EtaFile<Field>& file) {
  size_t n = file.get_breadth();

  for (const auto& entry : file) {
    auto matrix = CSCMatrix<Field>(n);

    for (size_t j = 0; j < n; ++j) {
      if (j != entry.index) {
        matrix.add_column();
        matrix.push_to_last_column(j, 1);
      } else {
        matrix.add_column(entry.values);
      }
    }

    os << matrix << "\n";
  }

  return os;
}
