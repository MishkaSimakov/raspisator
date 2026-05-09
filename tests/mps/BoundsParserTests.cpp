#include <gtest/gtest.h>

#include "mps/DataRecordTokenizer.h"
#include "mps/Types.h"
#include "mps/sections/BoundsParser.h"

using namespace mps;

TEST(BoundsParserTests, Simple) {
  MPSParsingState<double> state;
  BoundsParser<double> parser;

  state.add_col("x1");

  const auto string = " UP BOUND     x1                  40";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());
  parser.parse(record, state);

  ASSERT_TRUE(state.cols[0].upper_specified);
  ASSERT_FALSE(state.cols[0].lower_specified);
  ASSERT_EQ(state.cols[0].bound, (Bound<double>{0, 40}));
}

TEST(BoundsParserTests, BothLoAndUp) {
  MPSParsingState<double> state;
  BoundsParser<double> parser;

  state.add_col("JCH3TGBE");

  {
    const auto string = " LO BOUND     JCH3TGBE           10.";

    const auto record =
        DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());
    parser.parse(record, state);
  }

  {
    const auto string = " UP BOUND     JCH3TGBE           37.";

    const auto record =
        DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());
    parser.parse(record, state);
  }

  ASSERT_TRUE(state.cols[0].lower_specified);
  ASSERT_TRUE(state.cols[0].upper_specified);
  ASSERT_EQ(state.cols[0].bound, (Bound<double>{10, 37}));
}

TEST(BoundsParserTests, BothFxAndUp) {
  MPSParsingState<double> state;
  BoundsParser<double> parser;

  state.add_col("JCH3TGBE");

  {
    const auto string = " FX BOUND     JCH3TGBE           10.";

    const auto record =
        DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());
    parser.parse(record, state);
  }

  ASSERT_TRUE(state.cols[0].lower_specified);
  ASSERT_TRUE(state.cols[0].upper_specified);

  {
    const auto string = " UP BOUND     JCH3TGBE           37.";

    const auto record =
        DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());
    ASSERT_ANY_THROW({ parser.parse(record, state); });
  }
}

TEST(BoundsParserTests, BothLiAndUp) {
  MPSParsingState<double> state;
  BoundsParser<double> parser;

  state.add_col("JCH3TGBE");

  {
    const auto string = " LI BOUND     JCH3TGBE           10.";

    const auto record =
        DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());
    parser.parse(record, state);
  }

  {
    const auto string = " UP BOUND     JCH3TGBE           37.";

    const auto record =
        DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());
    parser.parse(record, state);
  }

  ASSERT_TRUE(state.cols[0].lower_specified);
  ASSERT_TRUE(state.cols[0].upper_specified);
  ASSERT_TRUE(state.cols[0].is_integer);
  ASSERT_EQ(state.cols[0].bound, (Bound<double>{10, 37}));
}

TEST(BoundsParserTests, DefaultBound) {
  MPSParsingState<double> state;
  BoundsParser<double> parser;

  state.add_col("JCH3TGBE");

  ASSERT_FALSE(state.cols[0].lower_specified);
  ASSERT_FALSE(state.cols[0].upper_specified);
  ASSERT_FALSE(state.cols[0].is_integer);
  ASSERT_EQ(state.cols[0].bound, (Bound<double>{0, std::nullopt}));
}

// specified in:
// https://www.ibm.com/docs/en/icos/22.1.0?topic=standard-records-in-mps-format
TEST(BoundsParserTests, DefaultBoundSwitch) {
  MPSParsingState<double> state;
  BoundsParser<double> parser;

  state.add_col("JCH3TGBE");

  const auto string = " UP BOUND     JCH3TGBE           -37.";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());
  parser.parse(record, state);

  ASSERT_FALSE(state.cols[0].lower_specified);
  ASSERT_TRUE(state.cols[0].upper_specified);
  ASSERT_FALSE(state.cols[0].is_integer);
  ASSERT_EQ(state.cols[0].bound, (Bound<double>{std::nullopt, -37}));
}
