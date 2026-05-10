#include <gtest/gtest.h>

#include <sstream>

#include "mps/MPSParser.h"

using namespace mps;

TEST(MPSParserTests, Simple) {
  std::string mps =
      "NAME test\n"
      "ROWS\n"
      " N obj\n"
      "COLUMNS\n"
      "   x1 obj 1\n"
      "RHS\n"
      "    RHS1 obj 0\n"
      "ENDATA";

  std::stringstream ss(mps);
  const auto state = MPSParser<double>::parse(ss, Format::FREE);

  ASSERT_EQ(state.problem_name, "test");
  ASSERT_EQ(state.rhs_vector_name, "RHS1");

  ASSERT_EQ(state.rows.size(), 1);
  ASSERT_EQ(state.rows[0].type, RowSense::FREE);
  ASSERT_EQ(state.rows[0].name, "obj");
  ASSERT_EQ(state.rows[0].rhs, 0);
  ASSERT_EQ(state.rows[0].range, std::nullopt);

  ASSERT_EQ(state.cols.size(), 1);
  ASSERT_EQ(state.cols[0].name, "x1");
  ASSERT_EQ(state.cols[0].bound, (Bound<double>{std::nullopt, std::nullopt}));
  ASSERT_EQ(state.cols[0].values.at(0), 1);
}

TEST(MPSParserTests, MissingSection) {
  // Each section of the MPS file except the RANGES and BOUNDS sections is
  // mandatory
  constexpr std::array mandatory_sections = {
      "NAME hello", "ROWS", "COLUMNS", "RHS", "ENDATA",
  };

  for (size_t i = 0; i < mandatory_sections.size(); ++i) {
    // omit i-th section
    std::string mps;
    for (size_t j = 0; j < mandatory_sections.size(); ++j) {
      if (j != i) {
        mps += mandatory_sections[j];
        mps += "\n";
      }
    }

    std::stringstream ss(mps);
    ASSERT_ANY_THROW({ MPSParser<double>::parse(ss, Format::FREE); });
  }
}

TEST(MPSParserTests, DuplicatedSection) {
  std::string mps =
      "NAME hello\n"
      "ROWS\n"
      "COLUMNS\n"
      "ROWS\n"
      "ENDATA";

  std::stringstream ss(mps);
  ASSERT_ANY_THROW({ MPSParser<double>::parse(ss, Format::FREE); });
}

TEST(MPSParserTests, MissingEndata) {
  std::string mps =
      "NAME test\n"
      "ROWS\n"
      "COLUMNS\n"
      "RHS\n";

  std::stringstream ss(mps);
  ASSERT_ANY_THROW({ MPSParser<double>::parse(ss, Format::FREE); });
}

// Duplicate NAME shouldn't cause error, because some SIF files contain
// duplicated NAME section (SCSD6.SIF)
TEST(MPSParserTests, DuplicateName) {
  std::string mps =
      "NAME test1\n"
      "NAME test2\n"
      "ROWS\n"
      "COLUMNS\n"
      "RHS\n"
      "ENDATA";

  std::stringstream ss(mps);
  auto state = MPSParser<double>::parse(ss, Format::FREE);

  ASSERT_EQ(state.problem_name, "test2");
}

// From QAP8.SIF
TEST(MPSParserTests, FixedFormatShortRow) {
  std::string mps =
      "NAME test\n"
      "ROWS\n"
      "  N NOBJ\n"
      "COLUMNS\n"
      "    Y001A001  NOBJ            10.0\n"
      "RHS\n"
      "ENDATA";

  std::stringstream ss(mps);
  auto state = MPSParser<double>::parse(ss, Format::FIXED);

  ASSERT_EQ(state.rows[0].name, "NOBJ    ");
  ASSERT_EQ(state.cols[0].values.size(), 1);
}
