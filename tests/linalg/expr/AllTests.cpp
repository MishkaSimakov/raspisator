#include <gtest/gtest.h>

#include "linalg/Vector.h"
#include "linalg/expr/All.h"
#include "linalg/expr/TransposedExpr.h"

using namespace linalg;

TEST(AllTests, WrapsInOwningView) {
  Vector<int> vec = {1, 2, 3};

  auto view = detail::all(std::move(vec));

  static_assert(std::same_as<decltype(view), detail::OwningView<Vector<int>>>);
}

TEST(AllTests, WrapsInRefView) {
  Vector<int> vec = {1, 2, 3};

  auto view = detail::all(vec);

  static_assert(std::same_as<decltype(view),
                             linalg::detail::RefView<linalg::Vector<int>>>);
}

TEST(AllTests, WrapsConstInRefView) {
  const Vector<int> vec = {1, 2, 3};

  auto view = detail::all(vec);

  static_assert(std::same_as<decltype(view), detail::RefView<Vector<int>>>);
}

TEST(AllTests, DoesNotWrapView) {
  const Vector<int> vec = {1, 2, 3};

  auto view = detail::all(detail::TransposedExpr(vec));

  static_assert(
      std::same_as<decltype(view),
                   detail::TransposedExpr<detail::RefView<Vector<int>>>>);
}
