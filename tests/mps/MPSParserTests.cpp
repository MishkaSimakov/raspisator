#include <gtest/gtest.h>

#include <sstream>

#include "mps/MPSParser.h"

using namespace mps;

TEST(MPSParserTests, Simple) {
  std::string mps =
      "NAME test\n"
      "ROWS\n"
      " N obj\n"
      "COLUMNS\n"
      "   x1 obj 1\n"
      "RHS\n"
      "    RHS1 obj 0\n"
      "ENDATA";

  std::stringstream ss(mps);

  const auto state = MPSParser<double>::parse(ss, Format::FREE);

  ASSERT_EQ(state.problem_name, "test");
  ASSERT_EQ(state.rhs_vector_name, "RHS1");

  ASSERT_EQ(state.rows.size(), 1);
  ASSERT_EQ(state.rows[0].type, RowSense::FREE);
  ASSERT_EQ(state.rows[0].name, "obj");
  ASSERT_EQ(state.rows[0].rhs, 0);
  ASSERT_EQ(state.rows[0].range, std::nullopt);

  ASSERT_EQ(state.cols.size(), 1);
  ASSERT_EQ(state.cols[0].name, "x1");
  ASSERT_EQ(state.cols[0].bound, (Bound<double>{0, std::nullopt}));
  ASSERT_EQ(state.cols[0].values.at(0), 1);
}
