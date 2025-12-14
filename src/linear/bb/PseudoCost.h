#pragma once

#include <vector>

#include "BaseAccountant.h"
#include "Settings.h"
#include "linear/bb/NodeRelativeLocation.h"
#include "linear/model/MILP.h"
#include "linear/simplex/BoundedSimplexMethod.h"
#include "utils/Accumulators.h"

template <typename Field, typename Accountant = BaseAccountant<Field>>
class PseudoCostBranchAndBound {
  size_t total_nodes_count_;

  // Node is pushed into waiting list when feasible solution in it is found and
  // branching variable is calculated
  struct Node {
    size_t id;

    std::vector<Field> upper;
    std::vector<Field> lower;

    std::vector<VariableState> variables_states;
    size_t branch_variable;

    Field value;
    Field expected;

    Field branch_value;
    Field fractional;
  };

  simplex::BoundedSimplexMethod<Field> lp_solver_;

  std::vector<Node> waiting_;
  std::optional<FiniteMILPSolution<Field>> incumbent_;

  std::vector<VariableType> variables_;

  std::vector<ArithmeticMean<Field>> upper_pseudocosts_;
  std::vector<ArithmeticMean<Field>> lower_pseudocosts_;

  BranchAndBoundSettings<Field> settings_;

  // used to calculate average simplex iterations count
  ArithmeticMean<double> simplex_iterations_;

  Accountant accountant_;

  static bool is_integer(const Field& value) {
    auto floored = FieldTraits<Field>::floor(value);

    return !FieldTraits<Field>::is_nonzero(floored - value);
  }

  bool is_candidate_incumbent(const Matrix<Field>& point) const {
    for (size_t i = 0; i < variables_.size(); ++i) {
      if (variables_[i] == VariableType::INTEGER && !is_integer(point[i, 0])) {
        return false;
      }
    }

    return true;
  }

  bool try_prune(const FiniteLPSolution<Field>& solution, const Node& node) {
    if (is_candidate_incumbent(solution.point)) {
      if (!incumbent_ || incumbent_->value < solution.value) {
        incumbent_ = FiniteMILPSolution<Field>(solution.point, solution.value);
      }

      accountant_.pruned_by_integrality(node.id, solution);

      return true;
    }

    if (incumbent_ && !FieldTraits<Field>::is_strictly_positive(
                          solution.value - incumbent_->value)) {
      accountant_.pruned_by_bounds(node.id, solution);
      return true;
    }

    return false;
  }

  // double score_branching_variable(const Node& node, const Matrix<Field>&
  // point,
  //                                 size_t index) const {
  //   // most infeasible branching
  //   // branch on the variable, whose fractional value is closest to 0.5
  //   Field f = point[index, 0] - FieldTraits<Field>::floor(point[index, 0]);
  //   return FieldTraits<Field>::abs(0.5 - f);
  // }

  Field merge_score(Field left_score, Field right_score) {
    return (1 - settings_.score_factor) * std::min(left_score, right_score) +
           settings_.score_factor * std::max(left_score, right_score);
  }

  std::optional<Field> score_via_simplex(
      size_t index, const Node& node, const FiniteLPSolution<Field>& solution) {
    // full strong branching
    // run simplex method for gamma iterations for each subproblem,
    Field fractional = FieldTraits<Field>::fractional(solution.point[index, 0]);
    Field branch_value = FieldTraits<Field>::floor(solution.point[index, 0]);

    // calculate left node
    std::optional<Field> lower_score;
    auto tight_upper = node.upper;
    tight_upper[index] = branch_value;

    try {
      lp_solver_.setup_warm_start(node.variables_states);
      auto run_result = lp_solver_.dual(node.lower, tight_upper);

      if (run_result.is_feasible()) {
        lower_score =
            (std::get<FiniteLPSolution<Field>>(run_result.solution).value -
             solution.value) /
            fractional;
        lower_pseudocosts_[index].record(*lower_score);
      }

      accountant_.strong_branching_simplex_run(node.id, index, run_result);
    } catch (...) {
      lp_solver_.setup_warm_start(node.variables_states);

      std::ofstream os(std::format("simplex_core_dump_{}.h", 0));
      lp_solver_.dump_state(os, 0, node.lower, tight_upper);
      throw;
    }

    // calculate right node
    std::optional<Field> upper_score;
    auto tight_lower = node.lower;
    tight_lower[index] = branch_value + 1;

    try {
      lp_solver_.setup_warm_start(node.variables_states);
      auto run_result = lp_solver_.dual(tight_lower, node.upper);

      if (run_result.is_feasible()) {
        upper_score =
            (std::get<FiniteLPSolution<Field>>(run_result.solution).value -
             solution.value) /
            (1 - fractional);
        upper_pseudocosts_[index].record(*upper_score);
      }

      accountant_.strong_branching_simplex_run(node.id, index, run_result);
    } catch (...) {
      lp_solver_.setup_warm_start(node.variables_states);

      std::ofstream os(std::format("simplex_core_dump_{}.h", 0));
      lp_solver_.dump_state(os, 0, tight_lower, node.upper);
      throw;
    }

    if (!lower_score || !upper_score) {
      return std::nullopt;
    }

    return merge_score(*lower_score, *upper_score);
  }

  // size_t find_branching_variable(const Node& node, const Matrix<Field>&
  // point) {
  //   std::optional<std::pair<size_t, double>> min_score;
  //
  //   for (size_t i = 0; i < variables_.size(); ++i) {
  //     if (variables_[i] == VariableType::INTEGER && !is_integer(point[i, 0]))
  //     {
  //       double score = score_branching_variable(node, point, i);
  //
  //       if (!min_score || score < min_score->second) {
  //         min_score = {i, score};
  //       }
  //     }
  //   }
  //
  //   assert(min_score.has_value());
  //   return min_score->first;
  // }

  bool is_reliable(size_t variable_id) {
    return std::min(lower_pseudocosts_[variable_id].count(),
                    upper_pseudocosts_[variable_id].count()) >=
           settings_.reliability_parameter;
  }

  Field score_via_pseudocost(size_t index, const Matrix<Field>& point) {
    Field fractional = FieldTraits<Field>::fractional(point[index, 0]);

    Field lower_cost = fractional * lower_pseudocosts_[index].mean();
    Field upper_cost = (1 - fractional) * upper_pseudocosts_[index].mean();

    return merge_score(lower_cost, upper_cost);
  }

  size_t find_branching_variable(const Node& node,
                                 const FiniteLPSolution<Field>& solution) {
    ArgMinimum<Field> score;

    // calculate pseudo cost score for each candidate
    for (size_t i = 0; i < variables_.size(); ++i) {
      if (variables_[i] == VariableType::INTEGER &&
          !is_integer(solution.point[i, 0])) {
        std::optional<Field> current_score =
            is_reliable(i) ? score_via_pseudocost(i, solution.point)
                           : score_via_simplex(i, node, solution);

        score.record(i, current_score);
      }
    }

    if (!score.argmin()) {
      // return first violating variable
      for (size_t i = 0; i < variables_.size(); ++i) {
        if (variables_[i] == VariableType::INTEGER &&
            !is_integer(solution.point[i, 0])) {
          return i;
        }
      }
    }

    assert(score.argmin().has_value());
    return *score.argmin();
  }

  void log_bounds(const Node& node,
                  const FiniteLPSolution<Field>& current_solution) {
    std::string lb = "*";

    if (incumbent_) {
      std::stringstream ss;
      ss << incumbent_->value;
      lb = ss.str();
    }

    std::println("#{}; LB: {}; UB: {}", node.id, lb, current_solution.value);

    size_t sum = 0;
    for (auto& cost : lower_pseudocosts_) {
      sum += cost.count();
    }
    std::cout << "sum: " << sum << std::endl;;
  }

  void record_simplex_run(const Node& parent, NodeRelativeLocation location,
                          const SimplexResult<Field>& run_result) {
    assert(location != NodeRelativeLocation::ROOT);

    simplex_iterations_.record(run_result.iterations_count);

    // TODO: should try to introduce some cost for infeasible problems
    if (std::holds_alternative<NoFeasibleElements>(run_result.solution)) {
      return;
    }

    const auto& solution =
        std::get<FiniteLPSolution<Field>>(run_result.solution);

    Field delta = solution.value - parent.value;

    if (location == NodeRelativeLocation::LEFT_CHILD) {
      lower_pseudocosts_[parent.branch_variable].record(delta /
                                                        parent.fractional);
    } else {
      upper_pseudocosts_[parent.branch_variable].record(
          delta / (1 - parent.fractional));
    }
  }

  void try_push_to_waiting(Node node, const SimplexResult<Field>& run_result) {
    // prune
    if (std::holds_alternative<NoFeasibleElements>(run_result.solution)) {
      accountant_.pruned_by_infeasibility(node.id);
      return;
    }

    auto solution = std::get<FiniteLPSolution<Field>>(run_result.solution);

    if (try_prune(solution, node)) {
      return;
    }

    log_bounds(node, solution);

    // branch
    // TODO: expected value
    node.value = solution.value;
    node.variables_states = solution.variables;

    node.branch_variable = find_branching_variable(node, solution);
    node.branch_value =
        FieldTraits<Field>::floor(solution.point[node.branch_variable, 0]);
    node.fractional =
        solution.point[node.branch_variable, 0] - node.branch_value;

    accountant_.set_branching_variable(node.id, solution, node.branch_variable,
                                       node.branch_value);

    waiting_.push_back(std::move(node));
  }

  void calculate_node(const Node& parent, NodeRelativeLocation location) {
    assert(location != NodeRelativeLocation::ROOT);

    Node child;
    child.id = ++total_nodes_count_;

    accountant_.set_child(parent.id, child.id, location);

    // tighten bounds
    child.lower = parent.lower;
    child.upper = parent.upper;

    if (location == NodeRelativeLocation::LEFT_CHILD) {
      child.upper[parent.branch_variable] = parent.branch_value;
    } else {
      child.lower[parent.branch_variable] = parent.branch_value + 1;
    }

    // run simplex method for subproblem
    lp_solver_.setup_warm_start(parent.variables_states);

    auto run_result = lp_solver_.dual(child.lower, child.upper);

    accountant_.simplex_run(child.id, run_result);
    record_simplex_run(parent, location, run_result);

    try_push_to_waiting(std::move(child), run_result);
  }

  static Matrix<Field> perturbed_costs(
      const MILPProblem<Field>& problem,
      const BranchAndBoundSettings<Field>& settings) {
    if (settings.perturbation == PerturbationMode::DISABLED) {
      return problem.c;
    }

    auto [n, d] = problem.A.shape();
    Matrix<Field> result = problem.c;

    if (settings.perturbation == PerturbationMode::CONSTANT) {
      for (size_t i = 0; i < d; ++i) {
        if (!FieldTraits<Field>::is_nonzero(problem.c[0, i])) {
          result[0, i] = settings.perturbation_value;
        }
      }
    } else if (settings.perturbation ==
               PerturbationMode::FOR_INTEGER_SOLUTION) {
      std::vector<size_t> perturbed;

      for (size_t i = 0; i < d; ++i) {
        if (FieldTraits<Field>::is_nonzero(problem.c[0, i])) {
          continue;
        }

        if (FieldTraits<Field>::is_nonzero(problem.lower_bounds[i])) {
          continue;
        }

        if (!FieldTraits<Field>::is_nonzero(problem.lower_bounds[i] -
                                            problem.upper_bounds[i])) {
          continue;
        }

        result[0, i] = static_cast<Field>(1) / problem.upper_bounds[i];
        perturbed.push_back(i);
      }

      for (size_t i : perturbed) {
        result[0, i] /= static_cast<Field>(2 * perturbed.size());
      }
    }

    return result;
  }

 public:
  explicit PseudoCostBranchAndBound(MILPProblem<Field> problem,
                                    BranchAndBoundSettings<Field> settings)
      : total_nodes_count_(0),
        lp_solver_(CSCMatrix(problem.A), problem.b,
                   perturbed_costs(problem, settings)),
        variables_(problem.variables_),
        upper_pseudocosts_(problem.variables_.size()),
        lower_pseudocosts_(problem.variables_.size()),
        settings_(settings) {
    Node root;

    root.id = ++total_nodes_count_;
    root.lower = problem.lower_bounds;
    root.upper = problem.upper_bounds;

    auto basic_variables = linalg::get_row_basis(linalg::transposed(problem.A));
    lp_solver_.setup_warm_start(basic_variables);
    root.variables_states = lp_solver_.get_variables_states();

    accountant_.set_root(root.id);
    waiting_.push_back(std::move(root));
  }

  const Accountant& get_accountant() const { return accountant_; }

  MILPSolution<Field> solve() {
    // solve for root
    Node root = waiting_.back();
    waiting_.pop_back();

    auto run_result = lp_solver_.dual(root.lower, root.upper);
    try_push_to_waiting(std::move(root), run_result);

    // main branch and bound cycle
    while (!waiting_.empty()) {
      if (settings_.max_nodes && total_nodes_count_ > settings_.max_nodes) {
        return ReachedNodesLimit{};
      }

      Node parent = waiting_.back();
      waiting_.pop_back();

      calculate_node(parent, NodeRelativeLocation::LEFT_CHILD);
      calculate_node(parent, NodeRelativeLocation::RIGHT_CHILD);
    }

    if (incumbent_) {
      return *incumbent_;
    }

    return NoFiniteSolution{};
  }
};
