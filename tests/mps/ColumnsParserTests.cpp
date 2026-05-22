#include <gtest/gtest.h>

#include "mps/Format.h"
#include "mps/detail/Types.h"
#include "mps/detail/sections/ColumnsParser.h"

using namespace mps;
using namespace mps::detail;

TEST(ColumnsParserTests, SimpleTest) {
  MPSParsingState<double> state;
  ColumnsParser<double> parser;

  state.add_row(RowSense::EQUAL, "COST");
  state.add_row(RowSense::EQUAL, "LIM1");

  const auto string =
      "    XONE      COST                 1   LIM1                 2";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());
  parser.parse(record, state);

  ASSERT_EQ(state.cols.size(), 1);
  ASSERT_EQ(state.cols_map.size(), 1);

  ASSERT_EQ(state.cols[0].name, "XONE");

  std::map<size_t, double> expected = {
      {0, 1},
      {1, 2},
  };

  ASSERT_EQ(state.cols[0].values, expected);
}

TEST(ColumnsParserTests, UnknownRowTest) {
  MPSParsingState<double> state;
  ColumnsParser<double> parser;

  const auto string =
      "    XONE      COST                 1   LIM1                 1";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());

  ASSERT_ANY_THROW({ parser.parse(record, state); });
}

TEST(ColumnsParserTests, MissingCoefficient1) {
  MPSParsingState<double> state;
  ColumnsParser<double> parser;

  state.add_row(RowSense::EQUAL, "COST");

  const auto string = "    XONE      COST                ";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());

  ASSERT_ANY_THROW({ parser.parse(record, state); });
}

TEST(ColumnsParserTests, MissingCoefficient2) {
  MPSParsingState<double> state;
  ColumnsParser<double> parser;

  state.add_row(RowSense::EQUAL, "COST");
  state.add_row(RowSense::EQUAL, "LIM1");

  const auto string =
      "    XONE      COST                 1   LIM1                 ";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());

  ASSERT_ANY_THROW({ parser.parse(record, state); });
}

TEST(ColumnsParserTests, NonConsecutiveColumns) {
  MPSParsingState<double> state;
  ColumnsParser<double> parser;

  state.add_row(RowSense::EQUAL, "COST1");
  state.add_row(RowSense::EQUAL, "COST2");
  state.add_row(RowSense::EQUAL, "COST3");

  {
    const auto string = "    XONE COST1 1";

    const auto record =
        DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());

    parser.parse(record, state);
  }

  {
    const auto string = "    XTWO COST2 1";

    const auto record =
        DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());

    parser.parse(record, state);
  }

  {
    const auto string = "    XONE COST3 1";

    const auto record =
        DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());

    ASSERT_ANY_THROW({ parser.parse(record, state); });
  }
}

TEST(ColumnsParserTests, DuplicateRow) {
  MPSParsingState<double> state;
  ColumnsParser<double> parser;

  state.add_row(RowSense::EQUAL, "COST1");

  const auto string = "    XONE COST1 1 COST1 2";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());

  ASSERT_ANY_THROW({ parser.parse(record, state); });
}

TEST(ColumnsParserTests, GarbageAfterValue) {
  MPSParsingState<double> state;
  ColumnsParser<double> parser;

  state.add_row(RowSense::EQUAL, "COST1");

  const auto string = "    XONE COST1 1abc";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());

  ASSERT_ANY_THROW({ parser.parse(record, state); });
}

TEST(ColumnsParserTests, GarbageInsteadOfValue) {
  MPSParsingState<double> state;
  ColumnsParser<double> parser;

  state.add_row(RowSense::EQUAL, "COST1");

  const auto string = "    XONE COST1 abc123";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());

  ASSERT_ANY_THROW({ parser.parse(record, state); });
}
