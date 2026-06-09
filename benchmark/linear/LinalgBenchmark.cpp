#include <benchmark/benchmark.h>

#include <random>

#include "linalg/CSCMatrix.h"
#include "linalg/Matrix.h"
#include "linalg/Transpose.h"

using namespace linalg;

static size_t N = 1'000;

CSCMatrix<double> random_sparse(size_t N, size_t average_per_col) {
  std::default_random_engine random;
  std::bernoulli_distribution nonzero_distribution(
      static_cast<double>(average_per_col) / static_cast<double>(N));

  auto result = CSCMatrix<double>::zeros(N);

  for (size_t col = 0; col < N; ++col) {
    result.add_column();

    for (size_t row = 0; row < N; ++row) {
      if (nonzero_distribution(random)) {
        result.push_to_last_column(row, 1);
      }
    }
  }

  return result;
}

static void NewLinalgLibrary(benchmark::State& state) {
  auto matrix = random_sparse(N, 10);

  auto cost = Matrix<double>::zeros(N, 1);
  auto pi = Matrix<double>::zeros(N, 1);

  for (auto _ : state) {
    Matrix result = cost - transpose(matrix) * pi;

    benchmark::DoNotOptimize(result);
    benchmark::ClobberMemory();
  }
}

static void OldWay(benchmark::State& state) {
  auto matrix = random_sparse(N, 10);

  auto cost = Matrix<double>::zeros(N, 1);
  auto pi = Matrix<double>::zeros(N, 1);

  for (auto _ : state) {
    auto result = Matrix<double>::zeros(N, 1);

    for (size_t i = 0; i < N; ++i) {
      result[i, 0] = cost[i, 0];

      for (const auto& [row, value] : matrix.col_entries(i)) {
        result[i, 0] -= value * pi[row, 0];
      }
    }

    benchmark::DoNotOptimize(result);
    benchmark::ClobberMemory();
  }
}

BENCHMARK(NewLinalgLibrary);
BENCHMARK(OldWay);

BENCHMARK_MAIN();
