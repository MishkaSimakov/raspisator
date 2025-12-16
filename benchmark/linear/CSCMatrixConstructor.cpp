#include <benchmark/benchmark.h>

#include <random>

#include "linear/sparse/CSCMatrix.h"

static size_t N = 1'000;

Matrix<double> sparse_matrix_1(size_t N, size_t density_multiplier) {
  auto A = Matrix<double>::unity(N);

  std::mt19937 engine(0);
  std::uniform_int_distribution<size_t> uniform_dist(0, N - 1);

  // add O(N) non-zeros in "random" places
  for (size_t i = 0; i < density_multiplier * N; ++i) {
    size_t row = uniform_dist(engine);
    size_t col = uniform_dist(engine);

    A[row, col] = i + 1;
  }

  return A;
}

static void AllocateDuringConstruction(benchmark::State& state) {
  auto matrix = sparse_matrix_1(N, state.range(0));

  for (auto _ : state) {
    std::vector<double> data;
    std::vector<size_t> indices;
    std::vector<size_t> index_pointers(matrix.get_width() + 1);

    auto [n, d] = matrix.shape();

    index_pointers[0] = 0;
    size_t nonzero_cnt = 0;

    for (size_t col = 0; col < d; ++col) {
      for (size_t row = 0; row < n; ++row) {
        if (FieldTraits<double>::is_nonzero(matrix[row, col])) {
          data.push_back(matrix[row, col]);
          indices.push_back(row);
          ++nonzero_cnt;
        }
      }

      index_pointers[col + 1] = nonzero_cnt;
    }

    benchmark::DoNotOptimize(data);
    benchmark::DoNotOptimize(indices);
    benchmark::DoNotOptimize(index_pointers);

    benchmark::ClobberMemory();
  }
}

static void AllocateThenConstruct(benchmark::State& state) {
  auto matrix = sparse_matrix_1(N, state.range(0));

  for (auto _ : state) {
    benchmark::DoNotOptimize(matrix);

    size_t nonzero_cnt = 0;

    auto [n, d] = matrix.shape();
    for (size_t row = 0; row < n; ++row) {
      for (size_t col = 0; col < d; ++col) {
        if (FieldTraits<double>::is_nonzero(matrix[row, col])) {
          ++nonzero_cnt;
        }
      }
    }

    std::vector<double> data(nonzero_cnt);
    std::vector<size_t> indices(nonzero_cnt);
    std::vector<size_t> index_pointers(matrix.get_width() + 1);

    index_pointers[0] = 0;
    nonzero_cnt = 0;

    for (size_t col = 0; col < d; ++col) {
      for (size_t row = 0; row < n; ++row) {
        if (FieldTraits<double>::is_nonzero(matrix[row, col])) {
          data[nonzero_cnt] = matrix[row, col];
          indices[nonzero_cnt] = row;
          ++nonzero_cnt;
        }
      }

      index_pointers[col + 1] = nonzero_cnt;
    }

    benchmark::DoNotOptimize(data);
    benchmark::DoNotOptimize(indices);
    benchmark::DoNotOptimize(index_pointers);

    benchmark::ClobberMemory();
  }
}

BENCHMARK(AllocateDuringConstruction)->RangeMultiplier(2)->Range(1, N);
BENCHMARK(AllocateThenConstruct)->RangeMultiplier(2)->Range(1, N);

BENCHMARK_MAIN();
