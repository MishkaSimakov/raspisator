#include <gtest/gtest.h>

#include "mps/DataRecordTokenizer.h"

using namespace mps;

TEST(DataRecordTokenizerTests, FreeFormatWithoutEndSpaces) {
  std::string_view record = " UP BOUNDROW  DD1CLIND          .001";

  const auto result = DataRecordTokenizer::parse(record, Format::FREE, true);

  DataRecord expected;

  expected.fields[0] = record.substr(1, 2);
  expected.fields[1] = record.substr(4, 8);
  expected.fields[2] = record.substr(14, 8);
  expected.fields[3] = record.substr(32, 4);

  ASSERT_EQ(expected, result);
}

TEST(DataRecordTokenizerTests, FreeFormatWithEndSpaces) {
  std::string_view record = " UP YSBOUND   XI1918          22639.   ";

  const auto result = DataRecordTokenizer::parse(record, Format::FREE, true);

  DataRecord expected;

  expected.fields[0] = record.substr(1, 2);
  expected.fields[1] = record.substr(4, 7);
  expected.fields[2] = record.substr(14, 6);
  expected.fields[3] = record.substr(30, 6);

  ASSERT_EQ(expected, result);
}

TEST(DataRecordTokenizerTests, FreeFormatWithoutField1) {
  std::string_view record =
      "    JAN.71.   HLPSN3       63.991989   HLPSN4        428.1499   ";

  const auto result = DataRecordTokenizer::parse(record, Format::FREE, false);

  DataRecord expected;

  expected.fields[1] = record.substr(4, 7);
  expected.fields[2] = record.substr(14, 6);
  expected.fields[3] = record.substr(27, 9);
  expected.fields[4] = record.substr(39, 6);
  expected.fields[5] = record.substr(53, 8);

  ASSERT_EQ(expected, result);
}

TEST(DataRecordTokenizerTests, Comments1) {
  std::string_view record = " L  MX        $ Magnesium Maximum      lbs";

  const auto result = DataRecordTokenizer::parse(record, Format::FREE, true);

  DataRecord expected;

  expected.fields[0] = record.substr(1, 1);
  expected.fields[1] = record.substr(4, 2);

  expected.comment = record.substr(15, 27);

  ASSERT_EQ(expected, result);
}

TEST(DataRecordTokenizerTests, FreeFormatComments2) {
  std::string_view record =
      "    VAR2      ROW1      4              $ ROW2      10.1";

  const auto result = DataRecordTokenizer::parse(record, Format::FREE, false);

  DataRecord expected;

  expected.fields[1] = record.substr(4, 4);
  expected.fields[2] = record.substr(14, 4);
  expected.fields[3] = record.substr(24, 1);

  expected.comment = record.substr(40, 15);

  ASSERT_EQ(expected, result);
}

TEST(DataRecordTokenizerTests, FixedFormatComments1) {
  std::string_view record = " L  MX        $ Magnesium Maximum      lbs";

  const auto result = DataRecordTokenizer::parse(record, Format::FIXED, true);

  DataRecord expected;

  expected.fields[0] = record.substr(1, 1);
  expected.fields[1] = record.substr(4, 8);

  expected.comment = record.substr(15, 27);

  ASSERT_EQ(expected, result);
}

TEST(DataRecordTokenizerTests, FixedFormatComments2) {
  std::string_view record =
      "    VAR2      ROW1      4              $ ROW2      10.1";

  const auto result = DataRecordTokenizer::parse(record, Format::FIXED, false);

  DataRecord expected;

  expected.fields[1] = record.substr(4, 8);
  expected.fields[2] = record.substr(14, 8);
  expected.fields[3] = record.substr(24, 12);

  expected.comment = record.substr(40, 15);

  ASSERT_EQ(expected, result);
}
