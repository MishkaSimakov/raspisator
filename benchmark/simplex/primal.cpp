#include <fstream>
#include <iostream>
#include <print>
#include <unordered_set>

#include "mps/MPS.h"
#include "presolve/Chain.h"
#include "presolve/passes/RemoveLinearlyDependentEqualities.h"
#include "presolve/passes/TransformToEqualities.h"
#include "simplex/Config.h"
#include "simplex/Simplex.h"
#include "simplex/init/primal/Phase1.h"
#include "utils/Paths.h"

#include "presolve/passes/Scaling.h"
#include "problem/StandardMILP.h"
#include "problem/mutations/RemoveRows.h"
#include "simplex/pricing/primal/SteepestEdge.h"

using Field = double;

// - QAP15 - stuck in degenerate iterations in phase 1
// - QAP12 - same as QAP15
//
// The remaining entries are the largest Netlib instances by rows*cols, which is
// the cost driver for PrimalSteepestEdge::recalculate_weights (a dense
// n x (d-n) block solve rebuilt every refactorization). They are skipped so the
// steepest-edge weight-error experiment stays tractable. QAP15/QAP12 also fall
// in this top-20 by size.
const std::set<std::string> skipped = {
    "QAP15",    "QAP12",    "STOCFOR3", "DFL001",  "FIT2P",
    "MAROS-R7", "80BAU3B",  "GREENBEA", "GREENBEB", "D2Q06C",
    "PILOT87",  "WOODW",    "TRUSS",    "BNL2",    "SHIP12L",
    "CYCLE",    "PILOT",    "STOCFOR2", "SCTAP3",  "SHIP08L"};

int main() {
  std::ofstream output(paths::log("benchmark_simplex_primal.csv"));

  if (!output) {
    throw std::runtime_error("Failed to open output file.");
  }

  std::println(
      output, "name,status,objective,phase1_iterations,phase2_iterations,time");

  auto problems_path = paths::resource("lp_problems");

  std::vector<std::filesystem::directory_entry> all_problems;
  std::ranges::copy(std::filesystem::directory_iterator{problems_path},
                    std::back_inserter(all_problems));

  for (size_t index = 0; index < all_problems.size(); ++index) {
    const auto path = all_problems[index].path();
    const auto name = path.filename().replace_extension("").string();

    if (path.extension() != ".SIF") {
      continue;
    }

    std::println("{:02}/{} {}", index + 1, all_problems.size(), name);

    if (skipped.contains(name)) {
      std::println("  SKIPPED");
      std::println(output, "{},{},{},{},{},{}", name, "SKIPPED", 0, 0, 0, 0);
      continue;
    }

    std::ifstream is(all_problems[index].path());
    if (!is) {
      std::println("  Failed to open problem file.");
      continue;
    }

    auto problem = mps::read<Field>(is, mps::Format::FIXED);
    std::println("  size: {} x {}", problem.matrix.rows(),
                 problem.matrix.cols());

    try {
      auto start = std::chrono::steady_clock::now();

      auto optimizer = presolve::Chain<Field>()
                           .add<presolve::Scaling<Field>>()
                           .add<presolve::TransformToEqualities<Field>>();

      problem::StandardLP standard_problem(optimizer.apply(problem));

      auto phase1 = simplex::primal_phase1(
          standard_problem,
          simplex::Config<Field>()
              .set_validate_input(true)
              .set_primal_pricing<simplex::PrimalSteepestEdge<Field>>()
              .set_max_iterations(100'000));

      if (!phase1) {
        std::println("  Failed to find primal feasible basis.");
        std::println(output, "{},{},{},{},{},{}", name, "PHASE1_ERROR", 0, 0, 0,
                     0);
        continue;
      }

      // standard_problem = problem::remove_rows(std::move(standard_problem),
                                              // phase1->redundant_rows);

      // auto solver = simplex::Simplex<Field>();
      //
      // solver.set_validate_input(true);
      // solver.set_problem(standard_problem);
      // solver.set_primal_pricing<simplex::PrimalSteepestEdge<Field>>();
      // solver.set_max_iterations(100'000);
      //
      // auto result = solver.primal(phase1->states);
      //
      // auto end = std::chrono::steady_clock::now();
      // std::println("  status: {}, time: {}, iterations: {}",
      //              to_string(result.status), end - start,
      //              result.iterations_count + phase1->iterations_count);
      //
      // Field objective = 0;
      //
      // if (result.status == simplex::Status::OPTIMAL) {
      //   // final objective value is negated because all problems in the
      //   // benchmark are minimization problems, but solver works with
      //   // maximization problems and cost coefficients are negated.
      //   objective =
      //       -linalg::dot(optimizer.inverse(solver.get_point()), problem.cost);
      //
      //   std::println("  objective: {}", objective);
      // }
      //
      // std::println(
      //     output, "{},{},{},{},{},{}", name, to_string(result.status),
      //     objective, phase1->iterations_count, result.iterations_count,
      //     std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
      //         .count());
    } catch (std::exception& error) {
      std::println(output, "{},{},{},{},{},{}", name, "EXCEPTION", 0, 0, 0, 0);
      std::println("  failed: {}", error.what());
    }

    output.flush();
  }

  return 0;
}
