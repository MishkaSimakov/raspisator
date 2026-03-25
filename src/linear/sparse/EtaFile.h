#pragma once

#include <type_traits>
#include <unordered_map>
#include <vector>

#include "CSCMatrix.h"
#include "linear/matrix/Matrix.h"
#include "utils/Accumulators.h"

template <typename T, bool is_const>
using const_if = std::conditional_t<is_const, const T, T>;

namespace linalg {

enum class EtaType { ROW, COLUMN };

template <typename Field>
class EtaFile {
  struct Entry {
    size_t begin;
    size_t index;
    EtaType type;

    bool is_removed;
  };

  template <bool is_const>
  struct EntryView {
    using ValueT = const_if<std::pair<size_t, Field>, is_const>;

    std::span<ValueT> values;

    const size_t index;
    const EtaType type;

    EntryView(std::span<ValueT> values, size_t index, EtaType type)
        : values(values), index(index), type(type) {}

    EntryView(EntryView<false> other)
      requires(is_const)
        : values(other.values), index(other.index), type(other.type) {}
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
          entry_->index, entry_->type);
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
    for (size_t col = 0; col < vector.get_width(); ++col) {
      if ((entry.type == EtaType::COLUMN) != transposed) {
        Field a = vector[entry.index, col];
        vector[entry.index, col] = 0;

        for (auto [row, value] : entry.values) {
          vector[row, col] += a * value;
        }
      } else {
        KahanSum<Field> dot;

        for (auto [row, value] : entry.values) {
          dot.add(value * vector[row, col]);
        }

        vector[entry.index, col] = dot.sum();
      }
    }

    return vector;
  }

 public:
  EtaFile() = default;

  void push_back(size_t pivot_index,
                 const std::vector<std::pair<size_t, Field>>& values,
                 EtaType type = EtaType::COLUMN) {
    entries_.push_back(Entry{
        .begin = values_.size(),
        .index = pivot_index,
        .type = type,
        .is_removed = false,
    });

    values_.insert(values_.end(), values.cbegin(), values.cend());
  }

  void push_back(size_t pivot_index, const Matrix<Field>& vector,
                 EtaType type = EtaType::COLUMN) {
    assert(vector.get_width() == 1);

    entries_.push_back(Entry{
        .begin = values_.size(),
        .index = pivot_index,
        .type = type,
        .is_removed = false,
    });

    for (size_t i = 0; i < vector.get_height(); ++i) {
      if (FieldTraits<Field>::is_nonzero(vector[i, 0])) {
        values_.emplace_back(i, vector[i, 0]);
      }
    }
  }

  Matrix<Field> apply(Matrix<Field> vector, EntryView<true> entry) const {
    return apply_impl(std::move(vector), entry, false);
  }

  Matrix<Field> apply_transposed(Matrix<Field> vector,
                                 EntryView<true> entry) const {
    return apply_impl(std::move(vector), entry, true);
  }

  CSCMatrix<Field> as_matrix(EntryView<true> entry) const {
    size_t n = get_breadth();
    auto result = CSCMatrix<Field>(n);

    if (entry.type == EtaType::COLUMN) {
      for (size_t j = 0; j < n; ++j) {
        if (j != entry.index) {
          result.add_column();
          result.push_to_last_column(j, 1);
        } else {
          result.add_column(entry.values);
        }
      }
    } else {
      std::unordered_map<size_t, Field> values;
      for (auto [col, value] : entry.values) {
        values.emplace(col, value);
      }

      for (size_t j = 0; j < n; ++j) {
        result.add_column();

        auto itr = values.find(j);

        if (itr == values.end() || j != entry.index) {
          result.push_to_last_column(j, 1);
        }

        if (itr != values.end()) {
          result.push_to_last_column(entry.index, itr->second);
        }
      }
    }

    return result;
  }

  void purge() {
    std::vector<std::pair<size_t, Field>> new_values;
    std::vector<Entry> new_entries;

    for (const auto& entry : *this) {
      new_entries.push_back(Entry{
          .begin = new_values.size(),
          .index = entry.index,
          .type = entry.type,
          .is_removed = false,
      });

      for (auto [index, value] : entry.values) {
        if (FieldTraits<Field>::is_nonzero(value)) {
          new_values.emplace_back(index, value);
        }
      }
    }

    values_ = std::move(new_values);
    entries_ = std::move(new_entries);
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

  size_t size() const { return entries_.size(); }

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
  for (const auto& entry : file) {
    os << file.as_matrix(entry) << "\n";
  }

  return os;
}
