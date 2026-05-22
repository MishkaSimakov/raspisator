#pragma once

#include "mps/detail/DataRecordTokenizer.h"

namespace mps::detail {

template <typename Field>
class SectionParser {
 public:
  virtual bool has_field_1() const = 0;

  virtual void parse(const DataRecord& record,
                     MPSParsingState<Field>& state) = 0;

  virtual void teardown() {}

  virtual ~SectionParser() = default;
};

}  // namespace mps::detail
