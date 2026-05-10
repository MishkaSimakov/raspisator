#include <gtest/gtest.h>

#include "mps/sections/ColumnsParser.h"

using namespace mps;

template <typename Field, typename Parser>
struct ParserTestWrapper {
  Parser parser;
  MPSParsingState<Field> state;
  Format format;

  explicit ParserTestWrapper(Format format) : format(format) {}

  friend ParserTestWrapper& operator<<(ParserTestWrapper& wrapper,
                                       std::string_view string) {
    const auto record = DataRecordTokenizer::parse(
        string, wrapper.format, wrapper.parser.has_field_1());

    wrapper.parser.parse(record, wrapper.state);

    return wrapper;
  }
};

TEST(IntegerVariablesExtensionTests, Simple) {
  ParserTestWrapper<double, ColumnsParser<double>> wrapper(Format::FREE);

  wrapper.state.add_row(RowSense::FREE, "obj");
  wrapper.state.add_row(RowSense::LESS_THAN, "c1");
  wrapper.state.add_row(RowSense::LESS_THAN, "c2");
  wrapper.state.add_row(RowSense::EQUAL, "c3");

  wrapper << " x1        obj                 -1   c1                  -1";
  wrapper << " MARK0000  'MARKER'                 'INTORG'";
  wrapper << " x4        obj                 -1   c1                  10";
  wrapper << " x4        c3                -3.5";
  wrapper << " MARK0001  'MARKER'                 'INTEND'";

  ASSERT_EQ(wrapper.state.cols[0].name, "x1");
  ASSERT_FALSE(wrapper.state.cols[0].is_integer);

  ASSERT_EQ(wrapper.state.cols[1].name, "x4");
  ASSERT_TRUE(wrapper.state.cols[1].is_integer);
}

TEST(IntegerVariablesExtensionTests, IntegerSectionNotClosed) {
  ParserTestWrapper<double, ColumnsParser<double>> wrapper(Format::FREE);

  wrapper.state.add_row(RowSense::FREE, "obj");
  wrapper.state.add_row(RowSense::LESS_THAN, "c1");
  wrapper.state.add_row(RowSense::LESS_THAN, "c2");
  wrapper.state.add_row(RowSense::EQUAL, "c3");

  wrapper << " x1        obj                 -1   c1                  -1";
  wrapper << " MARK0000  'MARKER'                 'INTORG'";
  wrapper << " x4        obj                 -1   c1                  10";
  wrapper << " x4        c3                -3.5";

  // should throw because integer section was never closed
  ASSERT_ANY_THROW({ wrapper.parser.teardown(); });
}

TEST(IntegerVariablesExtensionTests, ClosedBeforeOpened) {
  ParserTestWrapper<double, ColumnsParser<double>> wrapper(Format::FREE);

  wrapper.state.add_row(RowSense::FREE, "obj");
  wrapper.state.add_row(RowSense::FREE, "c1");

  wrapper << " x1        obj                 -1   c1                  -1";

  ASSERT_ANY_THROW(
      { wrapper << " MARK0000  'MARKER'                 'INTEND'"; });
}

TEST(IntegerVariablesExtensionTests, ClosedBeforeOpened2) {
  ParserTestWrapper<double, ColumnsParser<double>> wrapper(Format::FREE);

  wrapper.state.add_row(RowSense::FREE, "obj");
  wrapper.state.add_row(RowSense::FREE, "c1");

  wrapper << " x1        obj                 -1   c1                  -1";
  wrapper << " MARK0000  'MARKER'                 'INTORG'";
  wrapper << " MARK0001  'MARKER'                 'INTEND'";

  ASSERT_ANY_THROW(
      { wrapper << " MARK0002  'MARKER'                 'INTEND'"; });
}

TEST(IntegerVariablesExtensionTests, DoubleOpening) {
  ParserTestWrapper<double, ColumnsParser<double>> wrapper(Format::FREE);

  wrapper.state.add_row(RowSense::FREE, "obj");
  wrapper.state.add_row(RowSense::FREE, "c1");

  wrapper << " x1        obj                 -1   c1                  -1";
  wrapper << " MARK0000  'MARKER'                 'INTORG'";

  ASSERT_ANY_THROW(
      { wrapper << " MARK0001  'MARKER'                 'INTORG'"; });
}

TEST(IntegerVariablesExtensionTests, ColumnBothIntegerAndReal1) {
  ParserTestWrapper<double, ColumnsParser<double>> wrapper(Format::FREE);

  wrapper.state.add_row(RowSense::FREE, "obj");
  wrapper.state.add_row(RowSense::FREE, "c1");

  wrapper << " x1        obj                 -1";
  wrapper << " MARK0000  'MARKER'                 'INTORG'";

  ASSERT_ANY_THROW({ wrapper << "  x1        c1                  -1"; });
}

TEST(IntegerVariablesExtensionTests, ColumnBothIntegerAndReal2) {
  ParserTestWrapper<double, ColumnsParser<double>> wrapper(Format::FREE);

  wrapper.state.add_row(RowSense::FREE, "obj");
  wrapper.state.add_row(RowSense::FREE, "c1");

  wrapper << " MARK0000  'MARKER'                 'INTORG'";
  wrapper << "  x1        obj                  -1";
  wrapper << " MARK0000  'MARKER'                 'INTEND'";
  ASSERT_ANY_THROW({ wrapper << "  x1        c1                  -1"; });
}

TEST(IntegerVariablesExtensionTests, DuplicatedName1) {
  ParserTestWrapper<double, ColumnsParser<double>> wrapper(Format::FREE);

  wrapper.state.add_row(RowSense::FREE, "obj");
  wrapper.state.add_row(RowSense::FREE, "c1");

  wrapper << " x1        obj                 -1";
  wrapper << " x2        obj                 -1";
  ASSERT_ANY_THROW({ wrapper << " x1  'MARKER'                 'INTORG'"; });
}

TEST(IntegerVariablesExtensionTests, DuplicatedName2) {
  ParserTestWrapper<double, ColumnsParser<double>> wrapper(Format::FREE);

  wrapper.state.add_row(RowSense::FREE, "obj");
  wrapper.state.add_row(RowSense::FREE, "c1");

  wrapper << " x1        obj                 -1";
  wrapper << " x2        obj                 -1";
  wrapper << " MARK001  'MARKER'                 'INTORG'";
  wrapper << " MARK002  'MARKER'                 'INTEND'";
  wrapper << " x3        obj                 -1";
  ASSERT_ANY_THROW({ wrapper << " MARK002   obj                 2"; });
}

TEST(IntegerVariablesExtensionTests, EmptyIntegerSection) {
  ParserTestWrapper<double, ColumnsParser<double>> wrapper(Format::FREE);

  wrapper.state.add_row(RowSense::FREE, "obj");
  wrapper.state.add_row(RowSense::FREE, "c1");

  wrapper << " x1        obj                 -1";
  wrapper << " x2        obj                 -1";
  wrapper << " MARK001  'MARKER'                 'INTORG'";
  wrapper << " MARK002  'MARKER'                 'INTEND'";
  wrapper << " x3        obj                 -1";

  ASSERT_FALSE(wrapper.state.cols[0].is_integer);
  ASSERT_FALSE(wrapper.state.cols[1].is_integer);
  ASSERT_FALSE(wrapper.state.cols[2].is_integer);
}
