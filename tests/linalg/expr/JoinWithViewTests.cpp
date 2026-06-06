#include <gtest/gtest.h>

#include "linalg/expr/JoinWithView.h"

using linalg::detail::JoinWithView;

static_assert(
    std::ranges::range<JoinWithView<std::span<size_t>, std::span<size_t>>>);

TEST(JoinWithViewTests, SimpleTest) {
  std::vector<size_t> a = {1, 2, 3};
  std::vector<size_t> b = {4, 5, 6, 7};

  auto view = JoinWithView(a, b);

  std::vector result(view.begin(), view.end());
  std::vector<size_t> expected = {1, 2, 3, 4, 5, 6, 7};

  ASSERT_EQ(result, expected);
}

TEST(JoinWithViewTests, NestedView) {
  std::vector<size_t> a = {1, 2, 3};

  auto view = JoinWithView(
      a, a | std::views::transform([](size_t i) { return i + 3; }));

  std::vector result(view.begin(), view.end());
  std::vector<size_t> expected = {1, 2, 3, 4, 5, 6};

  ASSERT_EQ(result, expected);
}

// TEST(JoinWithViewTests, BidirectionalWalk) {
//   std::vector<size_t> a = {1, 2, 3};
//   std::vector<size_t> b = {4, 5, 6};
//
//   auto view = JoinWithView(a, b);
//
//   auto itr = std::ranges::begin(view);
//
//   ASSERT_EQ(*itr, 1);
//   ASSERT_EQ(*(++itr), 2);
//   ASSERT_EQ(*(++itr), 3);
//   ASSERT_EQ(*(++itr), 4);
//   ASSERT_EQ(*(++itr), 5);
//
//   ASSERT_EQ(*(--itr), 4);
//   ASSERT_EQ(*(--itr), 3);
//   ASSERT_EQ(*(--itr), 2);
//   ASSERT_EQ(*(--itr), 1);
//
//   ASSERT_EQ(*(++itr), 2);
//   ASSERT_EQ(*(++itr), 3);
//   ASSERT_EQ(*(++itr), 4);
//   ASSERT_EQ(*(++itr), 5);
//   ASSERT_EQ(*(++itr), 6);
//
//   ++itr;
//
//   ASSERT_EQ(itr, std::ranges::end(view));
// }

TEST(JoinWithViewTests, WithDropView) {
  std::vector<size_t> a = {1, 2, 3};
  std::set<size_t> b = {1, 2, 3, 4, 5, 6};

  auto view = JoinWithView(a, b | std::views::drop(3));

  std::vector result(view.begin(), view.end());
  std::vector<size_t> expected = {1, 2, 3, 4, 5, 6};

  ASSERT_EQ(result, expected);
}

TEST(JoinWithViewTests, OwningView) {
  std::vector<size_t> a = {1, 2, 3};
  std::set<size_t> b = {4, 5, 6};

  auto view = JoinWithView(a, std::move(b));

  auto other = std::move(view);

  std::vector result(other.begin(), other.end());
  std::vector<size_t> expected = {1, 2, 3, 4, 5, 6};

  ASSERT_EQ(result, expected);
}

TEST(JoinWithViewTests, IotaViewSecond) {
  std::vector<size_t> a = {1, 2, 3};
  auto b = std::views::iota(size_t{4});

  auto view = JoinWithView(a, b) | std::views::take(6);

  std::vector result( std::from_range, view);
  std::vector<size_t> expected = {1, 2, 3, 4, 5, 6};

  ASSERT_EQ(result, expected);
}
