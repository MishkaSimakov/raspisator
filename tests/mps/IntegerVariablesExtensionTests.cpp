#include <gtest/gtest.h>

#include "mps/sections/ColumnsParser.h"

using namespace mps;

template <typename Field>
void apply_parser(SectionParser<Field>& parser, MPSParsingState<Field>& state,
                  std::string string) {
  const auto record =
      DataRecordTokenizer::parse(string, Format::FREE, parser.has_field_1());
  parser.parse(record, state);
}

TEST(IntegerVariablesExtensionTests, Simple) {
  MPSParsingState<double> state;
  ColumnsParser<double> parser;

  state.add_row(RowSense::FREE, "obj");
  state.add_row(RowSense::LESS_THAN, "c1");
  state.add_row(RowSense::LESS_THAN, "c2");
  state.add_row(RowSense::EQUAL, "c3");

  apply_parser(parser, state,
               " x1        obj                 -1   c1                  -1");

  apply_parser(parser, state, " MARK0000  'MARKER'                 'INTORG'");

  apply_parser(parser, state,
               " x4        obj                 -1   c1                  10");

  apply_parser(parser, state, " x4        c3                -3.5");
  apply_parser(parser, state, " MARK0001  'MARKER'                 'INTEND'");

  ASSERT_EQ(state.cols[0].name, "x1");
  ASSERT_FALSE(state.cols[0].is_integer);

  ASSERT_EQ(state.cols[1].name, "x4");
  ASSERT_TRUE(state.cols[1].is_integer);
}
