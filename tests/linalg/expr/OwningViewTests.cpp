#include <gtest/gtest.h>

#include "linalg/MatrixSpy.h"
#include "linalg/expr/OwningView.h"

using namespace linalg;

TEST(OwningViewTests, CallsMoveConstructor) {
  MatrixSpy::reset_counters();

  MatrixSpy spy;
  auto view = detail::OwningView(std::move(spy));

  ASSERT_EQ(MatrixSpy::copy_constructor_calls, 0);
}

TEST(OwningViewTests, NonCopyable) {
  static_assert(!std::is_copy_assignable_v<detail::OwningView<MatrixSpy>>);
  static_assert(!std::is_copy_constructible_v<detail::OwningView<MatrixSpy>>);
}

TEST(OwningViewTests, Movable) {
  MatrixSpy spy;
  auto view = detail::OwningView(std::move(spy));

  // move construct
  MatrixSpy::reset_counters();
  auto other = std::move(view);

  ASSERT_EQ(MatrixSpy::copy_constructor_calls, 0);
  ASSERT_EQ(MatrixSpy::copy_assignment_calls, 0);

  // move assign
  MatrixSpy::reset_counters();
  view = std::move(other);

  ASSERT_EQ(MatrixSpy::copy_constructor_calls, 0);
  ASSERT_EQ(MatrixSpy::copy_assignment_calls, 0);
}
