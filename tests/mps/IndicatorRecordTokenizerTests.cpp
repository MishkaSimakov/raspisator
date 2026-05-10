#include <gtest/gtest.h>

#include "mps/IndicatorRecordTokenizer.h"

using namespace mps;

TEST(IndicatorRecordTokenizerTests, Simple) {
  const std::string line = "NAME hello world";

  const auto record = IndicatorRecordTokenizer::parse(line);

  ASSERT_EQ(record.type, SectionType::NAME);
  ASSERT_EQ(record.data, line.substr(5, 11));
}

TEST(IndicatorRecordTokenizerTests, ManySpaces) {
  const std::string line = "NAME      hello world  ";

  const auto record = IndicatorRecordTokenizer::parse(line);

  ASSERT_EQ(record.type, SectionType::NAME);
  ASSERT_EQ(record.data, line.substr(10, 13));
}

TEST(IndicatorRecordTokenizerTests, UnknownSectionName) {
  const std::string line = "ABRACADABRA";

  ASSERT_ANY_THROW({ IndicatorRecordTokenizer::parse(line); });
}

TEST(IndicatorRecordTokenizerTests, AllSectionNames) {
  const std::array sections = {
      SectionType::NAME,   SectionType::ROWS,   SectionType::COLUMNS,
      SectionType::RHS,    SectionType::RANGES, SectionType::BOUNDS,
      SectionType::ENDATA,
  };

  for (const SectionType type : sections) {
    const auto string = section_type_to_string(type);

    const auto record = IndicatorRecordTokenizer::parse(string);

    ASSERT_EQ(record.type, type);
    ASSERT_TRUE(record.data.empty());
  }
}
