#include <gtest/gtest.h>
#include <gmpxx.h>

TEST(GMPTest, ItWorks) {
  mpz_class x = 123;

  x += 321;

  ASSERT_EQ(x, 444);
}