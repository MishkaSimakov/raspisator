#include <gtest/gtest.h>

#include "mps/Types.h"
#include "mps/sections/ColumnsParser.h"

using namespace mps;

TEST(ColumnsParserTests, SimpleTest) {
  MPSParsingState<double> state;
  ColumnsParser<double> parser;

  state.add_row(RowSense::EQUAL, "COST");
  state.add_row(RowSense::EQUAL, "LIM1");

  const auto string =
      "    XONE      COST                 1   LIM1                 1";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());
  parser.parse(record, state);

  ASSERT_EQ(state.cols.size(), 1);
  ASSERT_EQ(state.cols_map.size(), 1);

  ASSERT_EQ(state.cols[0].name, "XONE");

  std::vector<std::pair<size_t, double>> expected = {
      {0, 1},
      {1, 1},
  };

  ASSERT_EQ(state.cols[0].values, expected);
}
