#include <gtest/gtest.h>

#include <functional>
#include <span>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "utils/MergeSortedRange.h"

namespace {

template <typename T>
using SparseVector = std::vector<std::pair<size_t, T>>;

template <typename T>
using Entry = std::tuple<size_t, T, T>;

// Collects the merged range into a flat vector so it can be compared against
// the expected result with a single ASSERT_EQ.
template <typename Range>
auto collect(const Range& range) {
  using value_type = typename std::decay_t<Range>::iterator::value_type;
  std::vector<value_type> result;
  for (const auto& entry : range) {
    result.push_back(entry);
  }
  return result;
}

}  // namespace

TEST(MergeSortedRangeTests, BothEmpty) {
  SparseVector<int> x;
  SparseVector<int> y;

  MergeSortedRange<int> merged{x, y};

  ASSERT_TRUE(merged.begin() == merged.end());
  ASSERT_TRUE(collect(merged).empty());
}

TEST(MergeSortedRangeTests, LeftEmpty) {
  SparseVector<int> x;
  SparseVector<int> y = {{0, 5}, {2, 7}, {4, 9}};

  MergeSortedRange<int> merged{x, y};

  std::vector<Entry<int>> expected = {{0, 0, 5}, {2, 0, 7}, {4, 0, 9}};
  ASSERT_EQ(collect(merged), expected);
}

TEST(MergeSortedRangeTests, RightEmpty) {
  SparseVector<int> x = {{0, 5}, {2, 7}, {4, 9}};
  SparseVector<int> y;

  MergeSortedRange<int> merged{x, y};

  std::vector<Entry<int>> expected = {{0, 5, 0}, {2, 7, 0}, {4, 9, 0}};
  ASSERT_EQ(collect(merged), expected);
}

TEST(MergeSortedRangeTests, DisjointKeys) {
  SparseVector<int> x = {{0, 10}, {2, 20}, {4, 30}};
  SparseVector<int> y = {{1, 11}, {3, 21}, {5, 31}};

  MergeSortedRange<int> merged{x, y};

  std::vector<Entry<int>> expected = {{0, 10, 0}, {1, 0, 11}, {2, 20, 0},
                                      {3, 0, 21}, {4, 30, 0}, {5, 0, 31}};
  ASSERT_EQ(collect(merged), expected);
}

TEST(MergeSortedRangeTests, IdenticalKeys) {
  SparseVector<int> x = {{0, 10}, {1, 20}, {2, 30}};
  SparseVector<int> y = {{0, 1}, {1, 2}, {2, 3}};

  MergeSortedRange<int> merged{x, y};

  std::vector<Entry<int>> expected = {{0, 10, 1}, {1, 20, 2}, {2, 30, 3}};
  ASSERT_EQ(collect(merged), expected);
}

TEST(MergeSortedRangeTests, PartialOverlap) {
  SparseVector<int> x = {{0, 10}, {2, 20}};
  SparseVector<int> y = {{1, 5}, {2, 7}};

  MergeSortedRange<int> merged{x, y};

  std::vector<Entry<int>> expected = {{0, 10, 0}, {1, 0, 5}, {2, 20, 7}};
  ASSERT_EQ(collect(merged), expected);
}

TEST(MergeSortedRangeTests, UnevenLengthsTrailingLeft) {
  SparseVector<int> x = {{1, 1}, {3, 3}, {5, 5}, {7, 7}};
  SparseVector<int> y = {{3, 30}};

  MergeSortedRange<int> merged{x, y};

  std::vector<Entry<int>> expected = {
      {1, 1, 0}, {3, 3, 30}, {5, 5, 0}, {7, 7, 0}};
  ASSERT_EQ(collect(merged), expected);
}

TEST(MergeSortedRangeTests, UnevenLengthsTrailingRight) {
  SparseVector<int> x = {{3, 30}};
  SparseVector<int> y = {{1, 1}, {3, 3}, {5, 5}, {7, 7}};

  MergeSortedRange<int> merged{x, y};

  std::vector<Entry<int>> expected = {
      {1, 0, 1}, {3, 30, 3}, {5, 0, 5}, {7, 0, 7}};
  ASSERT_EQ(collect(merged), expected);
}

TEST(MergeSortedRangeTests, SingleElementEach) {
  SparseVector<int> x = {{2, 42}};
  SparseVector<int> y = {{2, 99}};

  MergeSortedRange<int> merged{x, y};

  std::vector<Entry<int>> expected = {{2, 42, 99}};
  ASSERT_EQ(collect(merged), expected);
}

TEST(MergeSortedRangeTests, StructuredBindingsIteration) {
  SparseVector<int> x = {{0, 10}, {2, 20}};
  SparseVector<int> y = {{1, 5}, {2, 7}};

  MergeSortedRange<int> merged{x, y};

  std::vector<size_t> keys;
  int sum = 0;
  for (auto [key, xv, yv] : merged) {
    keys.push_back(key);
    sum += xv + yv;
  }

  ASSERT_EQ(keys, (std::vector<size_t>{0, 1, 2}));
  ASSERT_EQ(sum, 10 + 5 + 20 + 7);
}

TEST(MergeSortedRangeTests, DoubleValues) {
  SparseVector<double> x = {{0, 1.5}, {2, 2.5}};
  SparseVector<double> y = {{1, 0.25}, {2, 0.5}};

  MergeSortedRange<double> merged{x, y};

  std::vector<Entry<double>> expected = {
      {0, 1.5, 0.0}, {1, 0.0, 0.25}, {2, 2.5, 0.5}};
  ASSERT_EQ(collect(merged), expected);
}

// With std::greater the inputs must be sorted in strictly descending key order.
TEST(MergeSortedRangeTests, GreaterComparatorDescending) {
  SparseVector<int> x = {{4, 10}, {2, 20}, {0, 30}};
  SparseVector<int> y = {{4, 1}, {3, 2}, {1, 3}};

  MergeSortedRange<int, std::greater<>> merged{x, y, std::greater<>{}};

  std::vector<Entry<int>> expected = {{4, 10, 1}, {3, 0, 2}, {2, 20, 0},
                                      {1, 0, 3}, {0, 30, 0}};
  ASSERT_EQ(collect(merged), expected);
}

// A stateless user-defined comparator that orders keys in descending order,
// exercising the [[no_unique_address]] custom-comparator path.
TEST(MergeSortedRangeTests, CustomStatelessComparator) {
  struct DescendingByKey {
    bool operator()(size_t a, size_t b) const { return a > b; }
  };

  SparseVector<int> x = {{5, 50}, {3, 30}, {1, 10}};
  SparseVector<int> y = {{4, 40}, {3, 33}};

  MergeSortedRange<int, DescendingByKey> merged{x, y, DescendingByKey{}};

  std::vector<Entry<int>> expected = {
      {5, 50, 0}, {4, 0, 40}, {3, 30, 33}, {1, 10, 0}};
  ASSERT_EQ(collect(merged), expected);
}

// CTAD as used by callers such as RRQR: the value type and comparator are
// deduced from already-typed spans.
TEST(MergeSortedRangeTests, ClassTemplateArgumentDeduction) {
  SparseVector<int> x = {{0, 10}, {2, 20}};
  SparseVector<int> y = {{1, 5}, {2, 7}};

  std::span<const std::pair<size_t, int>> xs{x};
  std::span<const std::pair<size_t, int>> ys{y};

  MergeSortedRange merged{xs, ys, std::less<>{}};

  std::vector<Entry<int>> expected = {{0, 10, 0}, {1, 0, 5}, {2, 20, 7}};
  ASSERT_EQ(collect(merged), expected);
}
