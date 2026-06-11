#pragma once

#include "Instance.h"
#include "Tag.h"
#include "catalog/TextbookProblems.h"

namespace faker {

template <typename Field>
std::vector<Instance<Field>> catalog(Tag mask = Tag{0}) {
  std::vector<Instance<Field>> result;

  auto add = [&](std::vector<TaggedInstance<Field>> instances) {
    for (auto& [tag, instance] : instances) {
      if ((tag | mask) == tag) {
        result.push_back(std::move(instance));
      }
    }
  };

  add(detail::textbook<Field>());

  return result;
}

}  // namespace faker
