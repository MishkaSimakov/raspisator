#include <gtest/gtest.h>

#include "mps/Types.h"
#include "mps/sections/RowsParser.h"

using namespace mps;

TEST(RowsParserTests, SimpleTest) {
  MPSParsingState<double> state;
  RowsParser<double> parser;

  const auto string = " L  ROW1";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FIXED, parser.has_field_1());
  parser.parse(record, state);

  ASSERT_EQ(state.rows.size(), 1);

  ASSERT_EQ(state.rows[0].type, RowSense::LESS_THAN);
  ASSERT_STREQ(state.rows[0].name.c_str(), "ROW1");

  ASSERT_EQ(state.rows_map.at("ROW1"), 0);
}
