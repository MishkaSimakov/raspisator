#include <gtest/gtest.h>

#include "mps/Types.h"
#include "mps/sections/RangesParser.h"

using namespace mps;

TEST(RangesParserTests, Simple) {
  MPSParsingState<double> state;
  RangesParser<double> parser;

  state.add_row(RowSense::EQUAL, "c1");
  state.add_row(RowSense::EQUAL, "c2");

  const auto string =
      "    rhs       c1                  20   c2                  30";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());
  parser.parse(record, state);

  ASSERT_EQ(state.rows[0].range, 20);
  ASSERT_EQ(state.rows[1].range, 30);
}

TEST(RangesParserTests, DuplicatedRow) {
  MPSParsingState<double> state;
  RangesParser<double> parser;

  state.add_row(RowSense::EQUAL, "c1");
  state.add_row(RowSense::EQUAL, "c2");

  const auto string =
      "    rhs       c1                  20   c1                  30";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());

  ASSERT_ANY_THROW({ parser.parse(record, state); });
}

TEST(RangesParserTests, MultipleVectors) {
  MPSParsingState<double> state;
  RangesParser<double> parser;

  state.add_row(RowSense::EQUAL, "c1");
  state.add_row(RowSense::EQUAL, "c2");

  {
    const auto string = "    rhs       c1                  20";

    const auto record =
        DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());

    parser.parse(record, state);
  }

  ASSERT_EQ(state.rows[0].range, 20);
  ASSERT_EQ(state.rows[1].range, std::nullopt);

  {
    const auto string = "    rhs2      c2                  30";

    const auto record =
        DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());

    parser.parse(record, state);
  }

  ASSERT_EQ(state.rows[0].range, 20);
  ASSERT_EQ(state.rows[1].range, std::nullopt);
}

TEST(RangesParserTests, PartialRecord) {
  MPSParsingState<double> state;
  RangesParser<double> parser;

  state.add_row(RowSense::EQUAL, "c1");
  state.add_row(RowSense::EQUAL, "c2");

  const auto string = "    rhs       c1                  20";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());

  parser.parse(record, state);

  ASSERT_EQ(state.rows[0].range, 20);
  ASSERT_EQ(state.rows[1].range, std::nullopt);
}

TEST(RangesParserTests, UnknownRow) {
  MPSParsingState<double> state;
  RangesParser<double> parser;

  state.add_row(RowSense::EQUAL, "c1");
  state.add_row(RowSense::EQUAL, "c2");

  const auto string = "    rhs       c3                  20";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());

  ASSERT_ANY_THROW({ parser.parse(record, state); });
}
