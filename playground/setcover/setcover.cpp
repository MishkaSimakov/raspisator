#include "Reader.h"
#include "Types.h"
#include "linear/problem/MILPProblem.h"
#include "linear/problem/ToMatrices.h"
#include "linear/problem/optimization/TransformToEqualities.h"
#include "linear/simplex/BoundedSimplexMethod.h"

MILPProblem<double> to_milp(const setcover::Problem& problem) {
  MILPProblem<double> result;

  std::vector<Variable<double>> variables;
  for (size_t i = 0; i < problem.sets.size(); ++i) {
    variables.push_back(result.new_variable(std::format("x_{}", i),
                                            VariableType::INTEGER, 0, 1));
  }

  for (size_t j = 0; j < problem.elements_count; ++j) {
    Expression<double> constraint;

    for (size_t i = 0; i < problem.sets.size(); ++i) {
      if (problem.sets[i].elements.contains(j)) {
        constraint += variables[i];
      }
    }

    result.add_constraint(constraint >= Expression<double>{1});
  }

  Expression<double> objective;
  for (size_t i = 0; i < problem.sets.size(); ++i) {
    objective += static_cast<double>(problem.sets[i].cost) * variables[i];
  }
  result.set_objective(-objective);

  return result;
}

int main() {
  auto problem = setcover::read_problem(paths::resource("setcover/sc_10000_0"));

  auto milp_problem = to_milp(problem);

  auto optimizer = TransformToEqualities<double>();
  auto optimized = optimizer.apply(milp_problem);
  auto matrices = to_matrices(optimized);

  std::cout << matrices.A.get_height() << " x " << matrices.A.get_width()
            << std::endl;

  // run simplex algorithm
  auto simplex = simplex::BoundedSimplexMethod(CSCMatrix(matrices.A),
                                               matrices.b, matrices.c);

  auto result = simplex.dual(matrices.lower, matrices.upper);
}
