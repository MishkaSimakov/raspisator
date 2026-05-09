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

TEST(RowsParserTests, WrongSenseString) {
  MPSParsingState<double> state;
  RowsParser<double> parser;

  const auto string = " X  ROW1";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FIXED, parser.has_field_1());

  ASSERT_ANY_THROW({ parser.parse(record, state); });
}

TEST(RowsParserTests, EmptySenseString) {
  MPSParsingState<double> state;
  RowsParser<double> parser;

  const auto string = "    ROW1";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FIXED, parser.has_field_1());

  ASSERT_ANY_THROW({ parser.parse(record, state); });
}

TEST(RowsParserTests, DuplicateRow) {
  MPSParsingState<double> state;
  RowsParser<double> parser;

  {
    const auto string = " N  ROW1";
    const auto record =
        DataRecordTokenizer::parse(string, Format::FIXED, parser.has_field_1());

    parser.parse(record, state);
  }

  {
    const auto string = " L  ROW1";
    const auto record =
        DataRecordTokenizer::parse(string, Format::FIXED, parser.has_field_1());

    ASSERT_ANY_THROW({ parser.parse(record, state); });
  }
}

TEST(RowsParserTests, MultipleDifferentRows) {
  MPSParsingState<double> state;
  RowsParser<double> parser;

  for (size_t i = 0; i < 100; ++i) {
    const std::string string = std::format(" G  ROW{}", i);
    const auto record =
        DataRecordTokenizer::parse(string, Format::FIXED, parser.has_field_1());

    parser.parse(record, state);
  }

  ASSERT_EQ(state.rows.size(), 100);
  ASSERT_EQ(state.rows_map.size(), 100);

  for (size_t i = 0; i < 100; ++i) {
    const std::string string = std::format("ROW{}", i);
    ASSERT_EQ(state.rows_map.at(string), i);
  }
}
