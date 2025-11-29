#include <benchmark/benchmark.h>

#include <random>

#include "linear/matrix/LU.h"
#include "linear/matrix/Matrix.h"
#include "linear/sparse/LU.h"

static size_t N = 1'000;

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
  auto matrix = sparse_matrix_1(N, state.range(0));

  for (auto _ : state) {
    auto [L, U, P] = linalg::sparse_lu(matrix);

    benchmark::DoNotOptimize(L);
    benchmark::DoNotOptimize(U);
  }
}

static void DenseMatrix(benchmark::State& state) {
  auto matrix = linalg::to_dense(sparse_matrix_1(N, state.range(0)));

  for (auto _ : state) {
    auto [L, U, P] = linalg::get_lup(matrix);

    benchmark::DoNotOptimize(L);
    benchmark::DoNotOptimize(U);
  }
}

BENCHMARK(SparseMatrix)->RangeMultiplier(2)->Range(1, N);
BENCHMARK(DenseMatrix)->RangeMultiplier(2)->Range(1, N);

BENCHMARK_MAIN();
