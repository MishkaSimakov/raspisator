#pragma once

#include <vector>

#include "faker/Instance.h"
#include "faker/Tag.h"
#include "linalg/Linalg.h"
#include "problem/StandardLP.h"

namespace faker::detail {

template <typename Field>
TaggedInstance<Field> textbook1() {
  problem::StandardLP<Field> problem(2, 4);

  problem.name = "textbook1";

  problem.matrix = CSCMatrix<Field>(Matrix<Field>{
      {1, -1, 1, 0},
      {2, 1, 0, 1},
  });

  problem.rhs = {1, 3};
  problem.cost = {2, 1, 1, -1};

  problem.var_bounds = {Bound<Field>{0, 10}, Bound<Field>{0, 10},
                        Bound<Field>{0, 10}, Bound<Field>{0, 10}};

  TaggedInstance<Field> result;

  result.instance.problem = problem::MILP(problem);
  result.instance.solution_type = SolutionType::FEASIBLE;
  result.instance.problem_type = ProblemType::StandardLP;
  result.instance.feasible_point = {0, 3, 4, 0};
  result.instance.has_linearly_dependent_rows = false;
  result.instance.optimal_objective = 7;

  result.tag = Tag::TINY | Tag::FEASIBLE | Tag::STANDARD_LP |
               Tag::KNOWN_OPTIMAL_OBJECTIVE | Tag::ALL_VARIABLES_BOUNDED;

  return result;
}

template <typename Field>
TaggedInstance<Field> textbook2() {
  problem::StandardLP<Field> problem(2, 4);

  problem.name = "textbook2";

  problem.matrix = CSCMatrix<Field>(Matrix<Field>{
      {1, 1, -1, 1},
      {1, 14, 10, -10},
  });

  problem.rhs = {2, 24};
  problem.cost = {1, 2, 3, -4};

  problem.var_bounds = {Bound<Field>{0, 10}, Bound<Field>{0, 10},
                        Bound<Field>{0, 10}, Bound<Field>{0, 10}};

  TaggedInstance<Field> result;

  result.instance.problem = problem::MILP(problem);
  result.instance.solution_type = SolutionType::FEASIBLE;
  result.instance.problem_type = ProblemType::StandardLP;
  result.instance.feasible_point = {4, 0, 2, 0};
  result.instance.has_linearly_dependent_rows = false;
  result.instance.optimal_objective = 10;

  result.tag = Tag::TINY | Tag::FEASIBLE | Tag::STANDARD_LP |
               Tag::KNOWN_OPTIMAL_OBJECTIVE | Tag::ALL_VARIABLES_BOUNDED;

  return result;
}

template <typename Field>
TaggedInstance<Field> textbook3() {
  problem::StandardLP<Field> problem(1, 2);

  problem.name = "textbook3";

  problem.matrix = CSCMatrix<Field>(Matrix<Field>{{1, 1}});
  problem.rhs = {1};
  problem.cost = {1, 2};

  problem.var_bounds = {Bound<Field>{0, 10}, Bound<Field>{0, 10}};

  TaggedInstance<Field> result;

  result.instance.problem = problem::MILP(problem);
  result.instance.solution_type = SolutionType::FEASIBLE;
  result.instance.problem_type = ProblemType::StandardLP;
  result.instance.feasible_point = {0, 1};
  result.instance.has_linearly_dependent_rows = false;
  result.instance.optimal_objective = 2;

  result.tag = Tag::TINY | Tag::FEASIBLE | Tag::STANDARD_LP |
               Tag::KNOWN_OPTIMAL_OBJECTIVE | Tag::ALL_VARIABLES_BOUNDED;

  return result;
}

template <typename Field>
TaggedInstance<Field> textbook4() {
  // y -> max
  // x = 1
  // x in (-5, 5)
  // y in (0, +inf)

  problem::StandardLP<Field> problem(1, 2);

  problem.name = "textbook4";

  problem.matrix = CSCMatrix<Field>(Matrix<Field>{{1, 0}});
  problem.rhs = {{1}};
  problem.cost = {{0, 1}};

  problem.var_bounds = {Bound<Field>{-5, 5}, Bound<Field>{0, std::nullopt}};

  TaggedInstance<Field> result;

  result.instance.problem = problem::MILP(problem);
  result.instance.solution_type = SolutionType::UNBOUNDED;
  result.instance.problem_type = ProblemType::StandardLP;
  result.instance.feasible_point = {1, 0};
  result.instance.has_linearly_dependent_rows = false;
  result.instance.optimal_objective = std::nullopt;

  result.tag = Tag::TINY | Tag::UNBOUNDED | Tag::STANDARD_LP;

  return result;
}

// Variation of textbook1 with different bounds
template <typename Field>
TaggedInstance<Field> textbook5() {
  problem::StandardLP<Field> problem(2, 4);

  problem.name = "textbook5";

  problem.matrix = CSCMatrix<Field>(Matrix<Field>{
      {1, -1, 1, 0},
      {2, 1, 0, 1},
  });

  problem.rhs = {1, 3};
  problem.cost = {2, 1, 1, -1};

  problem.var_bounds = {Bound<Field>{0, 1}, Bound<Field>{0, 3},
                        Bound<Field>{1, 10}, Bound<Field>{0, 10}};

  TaggedInstance<Field> result;

  result.instance.problem = problem::MILP(problem);
  result.instance.solution_type = SolutionType::FEASIBLE;
  result.instance.problem_type = ProblemType::StandardLP;
  result.instance.feasible_point = {0, 0, 1, 3};
  result.instance.has_linearly_dependent_rows = false;
  result.instance.optimal_objective = 7;

  result.tag = Tag::TINY | Tag::FEASIBLE | Tag::STANDARD_LP |
               Tag::KNOWN_OPTIMAL_OBJECTIVE | Tag::ALL_VARIABLES_BOUNDED;

  return result;
}

template <typename Field>
TaggedInstance<Field> textbook6() {
  problem::StandardLP<Field> problem(1, 6);

  problem.name = "textbook6";

  problem.matrix = CSCMatrix<Field>(Matrix<Field>{{1, 1, 1, 0, 0, 0}});

  problem.rhs = {1};
  problem.cost = {1, 0, 0, 1, 1, 1};

  problem.var_bounds = std::vector<Bound<Field>>(6, {0, 1});

  TaggedInstance<Field> result;

  result.instance.problem = problem::MILP(problem);
  result.instance.solution_type = SolutionType::FEASIBLE;
  result.instance.problem_type = ProblemType::StandardLP;
  result.instance.feasible_point = {1, 0, 0, 1, 1, 1};
  result.instance.has_linearly_dependent_rows = false;
  result.instance.optimal_objective = 4;

  result.tag = Tag::TINY | Tag::FEASIBLE | Tag::STANDARD_LP |
               Tag::KNOWN_OPTIMAL_OBJECTIVE | Tag::ALL_VARIABLES_BOUNDED;

  return result;
}

template <typename Field>
TaggedInstance<Field> textbook7() {
  problem::StandardLP<Field> problem(3, 6);

  problem.name = "textbook7";

  problem.matrix = CSCMatrix<Field>(Matrix<Field>{
      {2, -1, 2, 1, 0, 0},
      {2, -3, 1, 0, 1, 0},
      {-1, 1, -2, 0, 0, 1},
  });

  problem.rhs = {4, -5, -1};
  problem.cost = {1, -1, 1, 0, 0, 0};

  problem.var_bounds = std::vector<Bound<Field>>(6, {0, std::nullopt});

  TaggedInstance<Field> result;

  result.instance.problem = problem::MILP(problem);
  result.instance.solution_type = SolutionType::FEASIBLE;
  result.instance.problem_type = ProblemType::StandardLP;
  result.instance.feasible_point = {
      0, Field(14) / Field(5), Field(17) / Field(5), 0, 0, 3};
  result.instance.has_linearly_dependent_rows = false;
  result.instance.optimal_objective = Field(3) / Field(5);

  result.tag = Tag::TINY | Tag::FEASIBLE | Tag::STANDARD_LP |
               Tag::KNOWN_OPTIMAL_OBJECTIVE;

  return result;
}

template <typename Field>
TaggedInstance<Field> textbook8() {
  problem::StandardLP<Field> problem(3, 6);

  problem.name = "textbook8";

  problem.matrix = CSCMatrix<Field>(Matrix<Field>{
      {2, -1, -2, 1, 0, 0},
      {2, -3, -1, 0, 1, 0},
      {-1, 1, 1, 0, 0, 1},
  });

  problem.rhs = {4, -5, -1};
  problem.cost = {1, -1, 1, 0, 0, 0};

  problem.var_bounds = std::vector<Bound<Field>>(6, {0, std::nullopt});

  TaggedInstance<Field> result;

  result.instance.problem = problem::MILP(problem);
  result.instance.solution_type = SolutionType::INFEASIBLE;
  result.instance.problem_type = ProblemType::StandardLP;
  result.instance.feasible_point = std::nullopt;
  result.instance.has_linearly_dependent_rows = false;
  result.instance.optimal_objective = std::nullopt;

  result.tag = Tag::TINY | Tag::INFEASIBLE | Tag::STANDARD_LP;

  return result;
}

template <typename Field>
TaggedInstance<Field> textbook9() {
  problem::StandardLP<Field> problem(3, 5);

  problem.name = "textbook9";

  problem.matrix = CSCMatrix<Field>(Matrix<Field>{
      {-5, 6, -8, -2, -8},
      {6, -8, -8, -2, 7},
      {4, 0, 2, 4, 4},
  });

  problem.rhs = {61, 113, -36};
  problem.cost = {-7, -10, -7, 9, 9};

  problem.var_bounds = {Bound<Field>{-6, 8}, Bound<Field>{-7, 2},
                        Bound<Field>{-10, 0}, Bound<Field>{-4, 5},
                        Bound<Field>{-2, 5}};

  TaggedInstance<Field> result;

  result.instance.problem = problem::MILP(problem);
  result.instance.solution_type = SolutionType::FEASIBLE;
  result.instance.problem_type = ProblemType::StandardLP;
  result.instance.feasible_point = {-Field(17) / Field(7), -7, -10,
                                    -Field(2) / Field(7), -Field(9) / Field(7)};
  result.instance.has_linearly_dependent_rows = false;
  result.instance.optimal_objective = Field(1000) / Field(7);

  result.tag = Tag::TINY | Tag::FEASIBLE | Tag::STANDARD_LP |
               Tag::KNOWN_OPTIMAL_OBJECTIVE | Tag::ALL_VARIABLES_BOUNDED;

  return result;
}

template <typename Field>
TaggedInstance<Field> textbook10() {
  problem::StandardLP<Field> problem(1, 2);

  problem.name = "textbook10";

  // \max x
  // x + y = 5
  // x free
  // y \in [1, 2]
  problem.matrix = CSCMatrix<Field>(Matrix<Field>{
      {1, 1},
  });

  problem.rhs = {5};
  problem.cost = {1, 0};

  problem.var_bounds = {Bound<Field>{std::nullopt, std::nullopt},
                        Bound<Field>{1, 2}};

  TaggedInstance<Field> result;

  result.instance.problem = problem::MILP(problem);
  result.instance.solution_type = SolutionType::FEASIBLE;
  result.instance.problem_type = ProblemType::StandardLP;
  result.instance.feasible_point = {4, 1};
  result.instance.has_linearly_dependent_rows = false;
  result.instance.optimal_objective = 4;

  result.tag = Tag::TINY | Tag::FEASIBLE | Tag::STANDARD_LP |
               Tag::KNOWN_OPTIMAL_OBJECTIVE;

  return result;
}

template <typename Field>
TaggedInstance<Field> textbook11() {
  problem::StandardLP<Field> problem(1, 2);

  problem.name = "textbook11";

  // \max x
  // x + y = 5
  // x free
  // y free
  problem.matrix = CSCMatrix<Field>(Matrix<Field>{
      {1, 1},
  });

  problem.rhs = {5};
  problem.cost = {1, 0};

  problem.var_bounds = {Bound<Field>{std::nullopt, std::nullopt},
                        Bound<Field>{std::nullopt, std::nullopt}};

  TaggedInstance<Field> result;

  result.instance.problem = problem::MILP(problem);
  result.instance.solution_type = SolutionType::UNBOUNDED;
  result.instance.problem_type = ProblemType::StandardLP;
  result.instance.feasible_point = std::nullopt;
  result.instance.has_linearly_dependent_rows = false;
  result.instance.optimal_objective = std::nullopt;

  result.tag = Tag::TINY | Tag::UNBOUNDED | Tag::STANDARD_LP;

  return result;
}

template <typename Field>
std::vector<TaggedInstance<Field>> textbook() {
  std::vector<TaggedInstance<Field>> result;

  result.push_back(textbook1<Field>());
  result.push_back(textbook2<Field>());
  result.push_back(textbook3<Field>());
  result.push_back(textbook4<Field>());
  result.push_back(textbook5<Field>());
  result.push_back(textbook6<Field>());
  result.push_back(textbook7<Field>());
  result.push_back(textbook8<Field>());
  result.push_back(textbook9<Field>());
  result.push_back(textbook10<Field>());
  result.push_back(textbook11<Field>());

  return result;
}

}  // namespace faker::detail
