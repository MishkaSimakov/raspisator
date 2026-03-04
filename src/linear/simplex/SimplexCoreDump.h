#pragma once

#include <chrono>
#include <fstream>

#include "linear/matrix/Matrix.h"
#include "linear/model/LP.h"
#include "linear/simplex/Simplex.h"
#include "linear/sparse/CSCMatrix.h"
#include "utils/String.h"

namespace simplex {

template <typename Field>
class SimplexCoreDump {
  const CSCMatrix<Field>& A_;
  const Matrix<Field>& b_;
  const Matrix<Field>& c_;

  static size_t get_dump_id() {
    return std::chrono::system_clock::now().time_since_epoch() /
           std::chrono::milliseconds(1);
  }

 public:
  SimplexCoreDump(const CSCMatrix<Field>& a, const Matrix<Field>& b,
                  const Matrix<Field>& c)
      : A_(a), b_(b), c_(c) {}

  void dump_state(const IterationState<Field>& state) {
    size_t dump_id = get_dump_id();
    std::string dump_name = std::format("simplex_core_dump_{}.h", dump_id);

    std::ofstream os(dump_name);

    if (!os) {
      throw std::runtime_error("Failed to open file for simplex core dump");
    }

    os << "namespace SimplexDump_" << dump_id << " {\n";

    os << "Matrix<Field> A = {" << linalg::to_dense(A_) << "};\n";
    os << "Matrix<Field> b = {" << b_ << "};\n";
    os << "Matrix<Field> c = {" << c_ << "};\n";

    std::vector<std::string> string_bounds(state.bounds->size());
    for (size_t i = 0; i < state.bounds->size(); ++i) {
      std::string bound = "std::pair{";

      if ((*state.bounds)[i].lower) {
        bound += std::format("{}", *(*state.bounds)[i].lower);
      } else {
        bound += "std::nullopt";
      }

      bound += ",";

      if ((*state.bounds)[i].upper) {
        bound += std::format("{}", *(*state.bounds)[i].upper);
      } else {
        bound += "std::nullopt";
      }

      bound += "}";

      string_bounds[i] = bound;
    }

    os << "Bounds<Field> bounds = {" << str::join(string_bounds, ", ")
       << "};\n";

    os << "std::vector<VariableState> last_states = {";
    for (auto var : state.variables_states) {
      if (var == VariableState::BASIC) {
        os << "VariableState::BASIC, ";
      } else if (var == VariableState::AT_LOWER) {
        os << "VariableState::AT_LOWER, ";
      } else {
        os << "VariableState::AT_UPPER, ";
      }
    }
    os << "};\n";

    os << "};";

    os << "}\n";

    os << std::flush;

    std::println("Registered failed simplex run into {}.", dump_name);
  }
};

}  // namespace simplex
