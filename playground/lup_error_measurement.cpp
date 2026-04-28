#include <iostream>

#include "linear/problem/MPS.h"
#include "linear/problem/optimization/RemoveLinearlyDependentConstraints.h"
#include "linear/problem/optimization/Scaling.h"
#include "linear/problem/optimization/TransformToEqualities.h"
#include "linear/simplex/Simplex.h"
#include "linear/sparse/CSCMatrix.h"
#include "linear/sparse/LU.h"
#include "linear/sparse/Norms.h"

MILPProblem<double> get_problem(const std::filesystem::path& path) {
  auto reader = MPSReader<double>(MPSFieldsMode::FIXED_WIDTH);
  reader.read(path);

  auto problem = reader.get_canonical_representation();

  problem = Scaling<double>().apply(problem);
  problem = TransformToEqualities<double>().apply(problem);
  problem = RemoveLinearlyDependentConstraints<double>().apply(problem);

  return problem;
}

int main() {
  const auto problem =
      to_matrices(get_problem(paths::resource("lp_problems/AFIRO.SIF")));
  const auto dense = problem.A;
  const auto sparse = CSCMatrix(dense);

  const auto [n, d] = sparse.shape();

  std::println("{} x {}", n, d);

  std::println("norm = {}", linalg::norm(sparse));

  auto lupa = linalg::LUPA<double>(sparse);

  std::vector<size_t> columns =
      linalg::get_row_basis(linalg::transposed(dense));
  lupa.set_columns(columns);

  const auto dense_submatrix = dense.get_columns(columns);

  Matrix<double> x(n, 1, 1);
  const auto Ax = dense_submatrix * x;

  const auto x_lupa = lupa.solve_linear(Ax);

  const auto residue = x - x_lupa;

  std::println("delta = {}", linalg::norm(residue));
}
