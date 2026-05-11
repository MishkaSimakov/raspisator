#include <gtest/gtest.h>

#include "mps/detail/Types.h"
#include "mps/detail/sections/RHSParser.h"

using namespace mps;
using namespace mps::detail;

TEST(RHSParserTests, Simple) {
  MPSParsingState<double> state;
  RHSParser<double> parser;

  state.add_row(RowSense::EQUAL, "c1");
  state.add_row(RowSense::EQUAL, "c2");

  const auto string =
      "    rhs       c1                  20   c2                  30";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());
  parser.parse(record, state);

  ASSERT_EQ(state.rows[0].rhs, 20);
  ASSERT_EQ(state.rows[1].rhs, 30);
}

TEST(RHSParserTests, DuplicatedRow) {
  MPSParsingState<double> state;
  RHSParser<double> parser;

  state.add_row(RowSense::EQUAL, "c1");
  state.add_row(RowSense::EQUAL, "c2");

  const auto string =
      "    rhs       c1                  20   c1                  30";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());

  ASSERT_ANY_THROW({ parser.parse(record, state); });
}

TEST(RHSParserTests, MultipleVectors) {
  MPSParsingState<double> state;
  RHSParser<double> parser;

  state.add_row(RowSense::EQUAL, "c1");
  state.add_row(RowSense::EQUAL, "c2");

  {
    const auto string = "    rhs       c1                  20";

    const auto record =
        DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());

    parser.parse(record, state);
  }

  ASSERT_EQ(state.rows[0].rhs, 20);
  ASSERT_EQ(state.rows[1].rhs, std::nullopt);

  {
    const auto string = "    rhs2      c2                  30";

    const auto record =
        DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());

    parser.parse(record, state);
  }

  ASSERT_EQ(state.rows[0].rhs, 20);
  ASSERT_EQ(state.rows[1].rhs, std::nullopt);
}

TEST(RHSParserTests, PartialRecord) {
  MPSParsingState<double> state;
  RHSParser<double> parser;

  state.add_row(RowSense::EQUAL, "c1");
  state.add_row(RowSense::EQUAL, "c2");

  const auto string = "    rhs       c1                  20";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());

  parser.parse(record, state);

  ASSERT_EQ(state.rows[0].rhs, 20);
  ASSERT_EQ(state.rows[1].rhs, std::nullopt);
}

TEST(RHSParserTests, UnknownRow) {
  MPSParsingState<double> state;
  RHSParser<double> parser;

  state.add_row(RowSense::EQUAL, "c1");
  state.add_row(RowSense::EQUAL, "c2");

  const auto string = "    rhs       c3                  20";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());

  ASSERT_ANY_THROW({ parser.parse(record, state); });
}
