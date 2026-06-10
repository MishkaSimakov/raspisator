#pragma once

#include <vector>

#include "SparseEtaMatrixView.h"
#include "linalg/Matrix.h"
#include "linalg/Vector.h"
#include "linear/FieldTraits.h"

namespace linalg {

template <typename Field>
class EtaFile {
  struct Entry {
    size_t begin;
    size_t index;
    EtaMatrixType type;

    bool is_removed;
  };

  size_t matrix_size_;
  std::vector<std::pair<size_t, Field>> values_;
  std::vector<Entry> entries_;

  template <bool is_const>
  class iterator_base {
    using EntryT = maybe_const_t<Entry, is_const>;
    using ValueT = maybe_const_t<std::pair<size_t, Field>, is_const>;

    using ValuesT =
        maybe_const_t<std::vector<std::pair<size_t, Field>>, is_const>;

    size_t matrix_size_;
    EntryT* entry_;
    EntryT* end_;
    ValuesT* values_;

    iterator_base(size_t matrix_size, EntryT* entry, EntryT* end,
                  ValuesT& values)
        : matrix_size_(matrix_size),
          entry_(entry),
          end_(end),
          values_(&values) {}

   public:
    using iterator_category = std::bidirectional_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = SparseEtaMatrixView<Field, is_const>;
    using reference = SparseEtaMatrixView<Field, is_const>;

    iterator_base() : entry_(nullptr), end_(nullptr), values_(nullptr) {}

    iterator_base(iterator_base<false> itr)
      requires(is_const)
        : entry_(itr.entry_), end_(itr.end_), values_(itr.values_) {}

    SparseEtaMatrixView<Field, is_const> operator*() const {
      size_t end = entry_ + 1 == end_ ? values_->size() : (entry_ + 1)->begin;

      return SparseEtaMatrixView<Field, is_const>(
          matrix_size_, entry_->index, entry_->type,
          std::span{values_->begin() + entry_->begin, values_->begin() + end});
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

  explicit EtaFile(size_t matrix_size) : matrix_size_(matrix_size) {}

  void push_back(size_t pivot_index,
                 const std::vector<std::pair<size_t, Field>>& values,
                 EtaMatrixType type = EtaMatrixType::COLUMN) {
    if (pivot_index >= matrix_size_) {
      throw std::invalid_argument(std::format(""));
    }

    entries_.push_back(Entry{
        .begin = values_.size(),
        .index = pivot_index,
        .type = type,
        .is_removed = false,
    });

    values_.insert(values_.end(), values.cbegin(), values.cend());
  }

  void push_back(size_t pivot_index, const Vector<Field>& vector,
                 EtaMatrixType type = EtaMatrixType::COLUMN) {
    entries_.push_back(Entry{
        .begin = values_.size(),
        .index = pivot_index,
        .type = type,
        .is_removed = false,
    });

    // TODO: think about removing drop
    for (size_t i = 0; i < vector.size(); ++i) {
      if (FieldTraits<Field>::should_drop(vector[i])) {
        continue;
      }

      values_.emplace_back(i, vector[i]);
    }
  }

  void purge() {
    std::vector<std::pair<size_t, Field>> new_values;
    std::vector<Entry> new_entries;

    for (const auto& entry : *this) {
      new_entries.push_back(Entry{
          .begin = new_values.size(),
          .index = entry.pivot_index(),
          .type = entry.type(),
          .is_removed = false,
      });

      new_values.append_range(entry.pivot_entries());
    }

    values_ = std::move(new_values);
    entries_ = std::move(new_entries);
  }

  void clear() {
    values_.clear();
    entries_.clear();
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

    return {matrix_size_, current, end, values_};
  }
  const_iterator begin() const { return const_cast<EtaFile*>(this)->begin(); }
  const_iterator cbegin() const { return begin(); }

  iterator end() {
    Entry* end = entries_.data() + entries_.size();
    return {matrix_size_, end, end, values_};
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

  //
  size_t rows() const { return matrix_size_; }
  size_t cols() const { return matrix_size_; }
  std::pair<size_t, size_t> shape() const { return {rows(), cols()}; }
};

}  // namespace linalg

template <typename Field>
std::ostream& operator<<(std::ostream& os, const linalg::EtaFile<Field>& file) {
  for (const auto& entry : file) {
    os << file.as_matrix(entry) << "\n";
  }

  return os;
}
