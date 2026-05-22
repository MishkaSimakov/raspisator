#include <gtest/gtest.h>

#include "mps/detail/DataRecordTokenizer.h"
#include "mps/detail/Types.h"
#include "mps/detail/sections/BoundsParser.h"

using namespace mps;
using namespace mps::detail;

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
  ASSERT_EQ(state.cols[0].bound, (Bound<double>{std::nullopt, 40}));
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

TEST(BoundsParserTests, Default) {
  MPSParsingState<double> state;
  BoundsParser<double> parser;

  state.add_col("JCH3TGBE");

  ASSERT_FALSE(state.cols[0].lower_specified);
  ASSERT_FALSE(state.cols[0].upper_specified);
  ASSERT_FALSE(state.cols[0].is_integer);
  ASSERT_EQ(state.cols[0].bound, (Bound<double>{std::nullopt, std::nullopt}));
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

TEST(BoundsParserTests, NoBoundSwitchForPl) {
  MPSParsingState<double> state;
  BoundsParser<double> parser;

  state.add_col("JCH3TGBE");

  const auto string = " PL BOUND     JCH3TGBE ";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());
  parser.parse(record, state);

  ASSERT_FALSE(state.cols[0].lower_specified);
  ASSERT_TRUE(state.cols[0].upper_specified);
  ASSERT_FALSE(state.cols[0].is_integer);
  ASSERT_EQ(state.cols[0].bound, (Bound<double>{std::nullopt, std::nullopt}));
}

TEST(BoundsParserTests, UnknownVariable) {
  MPSParsingState<double> state;
  BoundsParser<double> parser;

  state.add_col("x");

  const auto string = " UP BOUND     JCH3TGBE           -37.";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());

  ASSERT_ANY_THROW({ parser.parse(record, state); });
}

TEST(BoundsParserTests, FR) {
  MPSParsingState<double> state;
  BoundsParser<double> parser;

  state.add_col("JCH3TGBE");

  const auto string = " FR BOUND     JCH3TGBE";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());
  parser.parse(record, state);

  ASSERT_TRUE(state.cols[0].lower_specified);
  ASSERT_TRUE(state.cols[0].upper_specified);
  ASSERT_FALSE(state.cols[0].is_integer);
  ASSERT_EQ(state.cols[0].bound, (Bound<double>{std::nullopt, std::nullopt}));
}

TEST(BoundsParserTests, BV) {
  MPSParsingState<double> state;
  BoundsParser<double> parser;

  state.add_col("JCH3TGBE");

  const auto string = " BV BOUND     JCH3TGBE";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());
  parser.parse(record, state);

  ASSERT_TRUE(state.cols[0].lower_specified);
  ASSERT_TRUE(state.cols[0].upper_specified);
  ASSERT_TRUE(state.cols[0].is_integer);
  ASSERT_EQ(state.cols[0].bound, (Bound<double>{0, 1}));
}

TEST(BoundsParserTests, PlAndMi) {
  MPSParsingState<double> state;
  BoundsParser<double> parser;

  state.add_col("JCH3TGBE");

  {
    const auto string = " PL BOUND     JCH3TGBE";

    const auto record =
        DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());
    parser.parse(record, state);
  }

  {
    const auto string = " MI BOUND     JCH3TGBE";

    const auto record =
        DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());
    parser.parse(record, state);
  }

  ASSERT_TRUE(state.cols[0].lower_specified);
  ASSERT_TRUE(state.cols[0].upper_specified);
  ASSERT_FALSE(state.cols[0].is_integer);
  ASSERT_EQ(state.cols[0].bound, (Bound<double>{std::nullopt, std::nullopt}));
}

TEST(BoundsParserTests, UI) {
  MPSParsingState<double> state;
  BoundsParser<double> parser;

  state.add_col("JCH3TGBE");

  const auto string = " UI BOUND     JCH3TGBE 42";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());
  parser.parse(record, state);

  ASSERT_FALSE(state.cols[0].lower_specified);
  ASSERT_TRUE(state.cols[0].upper_specified);
  ASSERT_TRUE(state.cols[0].is_integer);
  ASSERT_EQ(state.cols[0].bound, (Bound<double>{std::nullopt, 42}));
}

TEST(BoundsParserTests, MultipleBoundsVectors) {
  MPSParsingState<double> state;
  BoundsParser<double> parser;

  state.add_col("x1");
  state.add_col("x2");

  {
    const auto string = " UP BOUND1     x1 42";

    const auto record =
        DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());
    parser.parse(record, state);
  }

  {
    const auto string = " UP BOUND2     x2 43";

    const auto record =
        DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());
    parser.parse(record, state);
  }

  ASSERT_FALSE(state.cols[1].lower_specified);
  ASSERT_FALSE(state.cols[1].upper_specified);
  ASSERT_FALSE(state.cols[1].is_integer);
  ASSERT_EQ(state.cols[1].bound, (Bound<double>{std::nullopt, std::nullopt}));
}

TEST(BoundsParserTests, TooManyColumns1) {
  MPSParsingState<double> state;
  BoundsParser<double> parser;

  state.add_col("x1");
  state.add_col("x2");

  const auto string = " UP BOUND1 x1 42 x2";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());
  ASSERT_ANY_THROW({ parser.parse(record, state); });
}

TEST(BoundsParserTests, TooManyColumns2) {
  MPSParsingState<double> state;
  BoundsParser<double> parser;

  state.add_col("x1");
  state.add_col("x2");

  const auto string = " UP BOUND1 x1 42 x2 3";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());
  ASSERT_ANY_THROW({ parser.parse(record, state); });
}

TEST(BoundsParserTests, FreeBoundWithValue) {
  // Value may be specified for a free bound, but it is ignored
  MPSParsingState<double> state;
  BoundsParser<double> parser;

  state.add_col("x1");

  const auto string = " FR BOUND1 x1 42";

  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());
  parser.parse(record, state);

  ASSERT_TRUE(state.cols[0].lower_specified);
  ASSERT_TRUE(state.cols[0].upper_specified);
  ASSERT_FALSE(state.cols[0].is_integer);
  ASSERT_EQ(state.cols[0].bound, (Bound<double>{std::nullopt, std::nullopt}));
}
