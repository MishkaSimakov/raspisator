#include <benchmark/benchmark.h>

#include <random>

#include "linear/matrix/Matrix.h"
#include "linear/sparse/LU.h"

CSCMatrix<double> sparse_matrix_1(size_t N, size_t density_multiplier) {
  auto A = Matrix<double>::unity(N);

  std::mt19937 engine(0);
  std::uniform_int_distribution<size_t> uniform_dist(0, N - 1);

  // add O(N) non-zeros in "random" places
  for (size_t i = 0; i < density_multiplier * N; ++i) {
    size_t row = uniform_dist(engine);
    size_t col = uniform_dist(engine);

    A[row, col] = i + 1;
  }

  return CSCMatrix(A);
}

static void SparseMatrix(benchmark::State& state) {
  auto matrix = sparse_matrix_1(100, 5);

  for (auto _ : state) {
    auto [L, U] = linalg::sparse_lu(matrix);

    benchmark::DoNotOptimize(L);
    benchmark::DoNotOptimize(U);
  }
}

BENCHMARK(SparseMatrix);

BENCHMARK_MAIN();
