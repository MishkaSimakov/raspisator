#pragma once

#include "Instance.h"
#include "catalog/TextbookProblems.h"

namespace faker {

template <typename Field>
class CatalogQueryBuilder {
  std::optional<ProblemType> problem_type_;
  std::optional<SolutionType> solution_type_;
  std::optional<bool> linearly_dependent_rows_;
  std::optional<bool> all_variables_bounded_;
  std::optional<bool> know_optimal_objective_;

  static bool are_all_variables_bounded(const Instance<Field>& instance) {
    for (size_t i = 0; i < instance.problem.var_bounds.size(); ++i) {
      const auto bound = instance.problem.var_bounds[i];

      if (!bound.lower || !bound.upper) {
        return false;
      }
    }

    return true;
  }

  bool satisfy_filter(const Instance<Field>& instance) const {
    if (problem_type_ && instance.problem_type != *problem_type_) {
      return false;
    }

    if (solution_type_ && instance.solution_type != *solution_type_) {
      return false;
    }

    if (linearly_dependent_rows_ &&
        instance.has_linearly_dependent_rows != *linearly_dependent_rows_) {
      return false;
    }

    if (all_variables_bounded_ &&
        are_all_variables_bounded(instance) != *all_variables_bounded_) {
      return false;
    }

    if (know_optimal_objective_ && (instance.optimal_objective !=
                                    std::nullopt) == *know_optimal_objective_) {
      return false;
    }

    return true;
  }

 public:
  CatalogQueryBuilder& problem_type(std::optional<ProblemType> value) {
    problem_type_ = value;
    return *this;
  }

  CatalogQueryBuilder& solution_type(std::optional<SolutionType> value) {
    solution_type_ = value;
    return *this;
  }

  CatalogQueryBuilder& linearly_dependent_rows(std::optional<bool> value) {
    linearly_dependent_rows_ = value;
    return *this;
  }

  CatalogQueryBuilder& all_variables_bounded(std::optional<bool> value) {
    all_variables_bounded_ = value;
    return *this;
  }

  CatalogQueryBuilder& know_optimal_objective(std::optional<bool> value) {
    know_optimal_objective_ = value;
    return *this;
  }

  std::vector<Instance<Field>> all() const {
    std::vector<Instance<Field>> result;

    auto add = [&](std::vector<Instance<Field>> instances) {
      for (auto& instance : instances) {
        if (satisfy_filter(instance)) {
          result.push_back(std::move(instance));
        }
      }
    };

    add(detail::textbook<Field>());

    return result;
  }
};

template <typename Field>
CatalogQueryBuilder<Field> catalog() {
  return CatalogQueryBuilder<Field>();
}

}  // namespace faker
